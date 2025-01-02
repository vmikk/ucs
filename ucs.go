package main

import (
	"bufio"
	"compress/gzip"
	"flag"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/parquet-go/parquet-go"
	"github.com/parquet-go/parquet-go/compress/zstd"
)

// Add new type to store command options
type Options struct {
	inputFile  string
	outputFile string
	summary    bool
	mapOnly    bool
	splitSeqID bool
	removeDups bool
}

// Add new types
type UCRecord struct {
	RecordType             string
	ClusterNumber          uint32
	SeqLengthOrClusterSize uint32
	PercentIdentity        *float64
	Strand                 *byte
	UnusedField1           string
	UnusedField2           string
	Alignment              string
	Query                  string
	Target                 string
}

// Add new type for Parquet output
type ParquetRecord struct {
	RecordType             string   `parquet:"record_type"`
	ClusterNumber          uint32   `parquet:"cluster_number"`
	SeqLengthOrClusterSize uint32   `parquet:"seq_length_or_cluster_size"`
	PercentIdentity        *float64 `parquet:"percent_identity"`
	Strand                 *byte    `parquet:"strand"`
	UnusedField1           string   `parquet:"unused_field_1"`
	UnusedField2           string   `parquet:"unused_field_2"`
	Alignment              string   `parquet:"alignment"`
	Query                  string   `parquet:"query"`
	Target                 string   `parquet:"target"`
}

func main() {
	opts := parseFlags()

	input, err := openInputFile(opts.inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
		os.Exit(1)
	}
	defer input.Close()

	output, err := createOutputFile(opts.outputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
		os.Exit(1)
	}
	defer output.Close()

	if opts.summary {
		rowCount, uniqueQuerySequences, uniqueTargetSequences, err := summarizeUC(input, opts.inputFile, opts)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error processing file: %v\n", err)
			os.Exit(1)
		}
		err = writeSummary(output, rowCount, uniqueQuerySequences, uniqueTargetSequences)
	} else {
		records, err := processUCFile(input, opts.inputFile, opts)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error processing file: %v\n", err)
			os.Exit(1)
		}

		writer := bufio.NewWriter(output)
		err = writeUCRecords(writer, records, opts)
		writer.Flush()
	}

	if err != nil {
		fmt.Fprintf(os.Stderr, "Error writing output: %v\n", err)
		os.Exit(1)
	}
}

func parseFlags() Options {
	opts := Options{}

	// Define flag pairs with long and short forms
	flagPairs := []struct {
		long, short string
		value       interface{}
		usage       string
		def         interface{}
	}{
		{"input", "i", &opts.inputFile, "Input file (default: stdin)", "-"},
		{"output", "o", &opts.outputFile, "Output file (default: stdout)", "-"},
		{"summary", "s", &opts.summary, "Print summary statistics", false},
		{"map-only", "m", &opts.mapOnly, "Output only Query-OTU mapping", true},
		{"split-id", "", &opts.splitSeqID, "Split sequence IDs at semicolon", true},
		{"rm-dups", "", &opts.removeDups, "Remove duplicate entries", true},
	}

	// Register all flags
	for _, f := range flagPairs {
		switch v := f.value.(type) {
		case *string:
			flag.StringVar(v, f.long, f.def.(string), f.usage)
			if f.short != "" {
				flag.StringVar(v, f.short, f.def.(string), f.usage)
			}
		case *bool:
			flag.BoolVar(v, f.long, f.def.(bool), f.usage)
			if f.short != "" {
				flag.BoolVar(v, f.short, f.def.(bool), f.usage)
			}
		}
	}

	// Custom usage message
	flag.Usage = func() {
		fmt.Fprintf(flag.CommandLine.Output(), "Usage of %s:\n", os.Args[0])
		for _, f := range flagPairs {
			shortFlag := ""
			if f.short != "" {
				shortFlag = fmt.Sprintf("-%s, ", f.short)
			}
			fmt.Fprintf(flag.CommandLine.Output(), "  %s--%s\t%s\n", shortFlag, f.long, f.usage)
		}
	}

	flag.Parse()
	return opts
}

func openInputFile(fileName string) (*os.File, error) {
	if fileName == "-" {
		return os.Stdin, nil
	}
	return os.Open(fileName)
}

func createOutputFile(fileName string) (*os.File, error) {
	if fileName == "-" {
		return os.Stdout, nil
	}
	return os.Create(fileName)
}

func processUCFile(input *os.File, inputFileName string, opts Options) ([]UCRecord, error) {
	scanner, err := createScanner(input, inputFileName)
	if err != nil {
		return nil, err
	}

	var records []UCRecord
	seenPairs := make(map[string]struct{}) // For duplicate detection
	duplicateCount := 0                    // Track number of duplicates

	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.Split(line, "\t")
		if len(fields) < 10 {
			continue
		}

		// Skip C records - they are redundant with S records
		// H & S records go first in the file, while C records are added to the very end
		if fields[0] == "C" {
			continue
		}

		// Process labels directly
		queryLabel := splitSeqID(fields[8], opts.splitSeqID)
		targetLabel := splitSeqID(fields[9], opts.splitSeqID)

		record := UCRecord{
			RecordType:             fields[0],
			ClusterNumber:          parseIntOrZero(fields[1]),
			SeqLengthOrClusterSize: parseIntOrZero(fields[2]),
			PercentIdentity:        parseFloatOrZero(fields[3]),
			Strand:                 parseStrand(fields[4]),
			UnusedField1:           fields[5],
			UnusedField2:           fields[6],
			Alignment:              fields[7],
			Query:                  queryLabel,
			Target:                 targetLabel,
		}

		// Process OTU mapping based on record type
		switch record.RecordType {
		case "H":
			// Hit record - use target as OTU
			record.Target = targetLabel
		case "C":
			// Cluster centroid - use query as both query and target
			record.Target = queryLabel
		default:
			// Skip other record types (N records)
			continue
		}

		// Check for duplicates if removal is enabled
		if opts.removeDups {
			pairKey := record.Query + "\t" + record.Target
			if _, exists := seenPairs[pairKey]; exists {
				duplicateCount++
				continue
			}
			seenPairs[pairKey] = struct{}{}
		}

		records = append(records, record)
	}

	// Print duplicate count to stderr in red color if any duplicates were found
	if duplicateCount > 0 {
		fmt.Fprintf(os.Stderr, "\033[31mRemoved %d duplicate entries\033[0m\n", duplicateCount)
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("reading input: %w", err)
	}

	return records, nil
}

func writeUCRecords(writer *bufio.Writer, records []UCRecord, opts Options) error {
	// Check if output is Parquet format
	if strings.HasSuffix(opts.outputFile, ".parquet") {
		return writeParquetRecords(opts.outputFile, records, opts)
	}

	if opts.mapOnly {
		// Write header for Query-OTU mapping
		_, err := fmt.Fprintf(writer, "Query\tOTU\n")
		if err != nil {
			return err
		}

		// Write data
		for _, record := range records {
			_, err := fmt.Fprintf(writer, "%s\t%s\n", record.Query, record.Target)
			if err != nil {
				return err
			}
		}
	} else {
		// Write header for full table
		_, err := fmt.Fprintf(writer, "recordType\tclusterNumber\tseqLengthOrClusterSize\tpercentIdentity\tstrand\tunusedField1\tunusedField2\talignment\tquery\ttarget\n")
		if err != nil {
			return err
		}

		// Write data
		for _, record := range records {
			strandStr := "*"
			if record.Strand != nil {
				strandStr = string(*record.Strand)
			}

			identityStr := "*"
			if record.PercentIdentity != nil {
				identityStr = fmt.Sprintf("%.2f", *record.PercentIdentity)
			}

			_, err := fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
				record.RecordType,
				record.ClusterNumber,
				record.SeqLengthOrClusterSize,
				identityStr,
				strandStr,
				record.UnusedField1,
				record.UnusedField2,
				record.Alignment,
				record.Query,
				record.Target)
			if err != nil {
				return err
			}
		}
	}
	return nil
}

// Add new function for Parquet writing
func writeParquetRecords(outputFile string, records []UCRecord, opts Options) error {
	// Open the output file
	f, err := os.Create(outputFile)
	if err != nil {
		return fmt.Errorf("creating output file: %w", err)
	}
	defer f.Close()

	// Map-only mode (Query-OTU table)
	if opts.mapOnly {
		// Define simple record type for map-only mode
		type SimpleRecord struct {
			Query string `parquet:"query"`
			OTU   string `parquet:"otu"`
		}

		// Create writer with ZSTD compression
		writer := parquet.NewGenericWriter[SimpleRecord](f,
			parquet.Compression(&zstd.Codec{}),
		)
		defer writer.Close()

		// Convert to simplified records
		simpleRecords := make([]SimpleRecord, len(records))
		for i, record := range records {
			simpleRecords[i].Query = record.Query
			simpleRecords[i].OTU = record.Target
		}

		// Write records
		if _, err := writer.Write(simpleRecords); err != nil {
			return fmt.Errorf("writing parquet records: %w", err)
		}
	} else {
		// Create writer with ZSTD compression for full records
		writer := parquet.NewGenericWriter[ParquetRecord](f,
			parquet.Compression(&zstd.Codec{}),
		)
		defer writer.Close()

		// Convert to parquet records
		parquetRecords := make([]ParquetRecord, len(records))
		for i, record := range records {
			parquetRecords[i] = ParquetRecord{
				RecordType:             record.RecordType,
				ClusterNumber:          record.ClusterNumber,
				SeqLengthOrClusterSize: record.SeqLengthOrClusterSize,
				PercentIdentity:        record.PercentIdentity,
				Strand:                 record.Strand,
				UnusedField1:           record.UnusedField1,
				UnusedField2:           record.UnusedField2,
				Alignment:              record.Alignment,
				Query:                  record.Query,
				Target:                 record.Target,
			}
		}

		// Write records
		if _, err := writer.Write(parquetRecords); err != nil {
			return fmt.Errorf("writing parquet records: %w", err)
		}
	}

	return nil
}

// Helper functions
func splitSeqID(id string, split bool) string {
	if !split {
		return id
	}
	parts := strings.Split(id, ";")
	return parts[0]
}

func parseIntOrZero(s string) uint32 {
	val, err := strconv.ParseUint(s, 10, 32)
	if err != nil {
		return 0
	}
	return uint32(val)
}

func parseFloatOrZero(s string) *float64 {
	if s == "*" {
		return nil
	}
	val, err := strconv.ParseFloat(s, 64)
	if err != nil {
		return nil
	}
	return &val
}

func summarizeUC(input *os.File, inputFileName string, opts Options) (int, int, int, error) {
	scanner, err := createScanner(input, inputFileName)
	if err != nil {
		return 0, 0, 0, err
	}

	rowCount := 0
	querySequences := make(map[string]struct{})
	targetSequences := make(map[string]struct{})

	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.Split(line, "\t")
		if len(fields) < 10 {
			continue
		}

		// Skip S records as they are redundant
		if fields[0] == "S" {
			continue
		}

		rowCount++
		queryLabel := splitSeqID(fields[8], opts.splitSeqID)
		querySequences[queryLabel] = struct{}{}

		// For centroid records (C), use Query as Target
		if fields[0] == "C" {
			targetSequences[queryLabel] = struct{}{}
		} else if fields[9] != "*" {
			targetLabel := splitSeqID(fields[9], opts.splitSeqID)
			targetSequences[targetLabel] = struct{}{}
		}
	}

	if err := scanner.Err(); err != nil {
		return 0, 0, 0, fmt.Errorf("reading input: %w", err)
	}

	return rowCount, len(querySequences), len(targetSequences), nil
}

func writeSummary(output *os.File, rowCount, uniqueQuerySequences, uniqueTargetSequences int) error {
	writer := bufio.NewWriter(output)
	defer writer.Flush()

	_, err := fmt.Fprintf(writer, "Total records: %d\n", rowCount)
	if err != nil {
		return err
	}

	_, err = fmt.Fprintf(writer, "Unique query sequences: %d\n", uniqueQuerySequences)
	if err != nil {
		return err
	}

	_, err = fmt.Fprintf(writer, "Unique target sequences: %d\n", uniqueTargetSequences)
	if err != nil {
		return err
	}

	return nil
}

// Helper function to parse strand
func parseStrand(s string) *byte {
	if s == "*" || len(s) == 0 {
		return nil
	}
	b := s[0] // take first byte
	if b != '+' && b != '-' {
		return nil
	}
	return &b
}

// New helper function to avoid code duplication
func createScanner(input *os.File, inputFileName string) (*bufio.Scanner, error) {
	var scanner *bufio.Scanner

	// Check if input is gzipped, either by filename or content
	isGzipped := strings.HasSuffix(inputFileName, ".gz")
	if !isGzipped && inputFileName == "-" {
		// For stdin, peek at the first two bytes to check for gzip magic number
		reader := bufio.NewReader(input)
		magic, err := reader.Peek(2)
		if err == nil && len(magic) == 2 && magic[0] == 0x1f && magic[1] == 0x8b {
			isGzipped = true
		}
		scanner = bufio.NewScanner(reader)
	}

	if isGzipped {
		gzipReader, err := gzip.NewReader(input)
		if err != nil {
			return nil, fmt.Errorf("creating gzip reader: %w", err)
		}
		scanner = bufio.NewScanner(gzipReader)
	} else if scanner == nil {
		scanner = bufio.NewScanner(input)
	}

	return scanner, nil
}
