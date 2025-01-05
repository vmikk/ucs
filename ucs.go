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

// UC record type
type UCRecord struct {
	RecordType    string   // Field 0: Record type (C/S/H/N)
	ClusterNumber uint32   // Field 1: Cluster number
	Size          uint32   // Field 2: Sequence length/cluster size
	Identity      *float64 // Field 3: % identity with centroid
	Strand        *byte    // Field 4: Strand +/-
	Unused1       string   // Field 5: unused
	Unused2       string   // Field 6: unused
	CIGAR         string   // Field 7: CIGAR string
	Query         string   // Field 8: Query sequence ID
	Target        string   // Field 9: Target/centroid sequence ID
}

// Add new type for Parquet output
type ParquetRecord struct {
	RecordType    string   `parquet:"record_type"`
	ClusterNumber uint32   `parquet:"cluster_number"`
	Size          uint32   `parquet:"size"`
	Identity      *float64 `parquet:"identity"`
	Strand        *byte    `parquet:"strand"`
	Unused1       string   `parquet:"unused_1"`
	Unused2       string   `parquet:"unused_2"`
	CIGAR         string   `parquet:"cigar"`
	Query         string   `parquet:"query"`
	Target        string   `parquet:"target"`
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
		rowCount, uniqueQuerySequences, uniqueTargetSequences, multiMappedQueries, err := summarizeUC(input, opts.inputFile, opts)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error processing file: %v\n", err)
			os.Exit(1)
		}
		err = writeSummary(output, rowCount, uniqueQuerySequences, uniqueTargetSequences, multiMappedQueries)
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
		{"split-id", "", &opts.splitSeqID, "Split sequence IDs at semicolon (default: true)", true},
		{"rm-dups", "", &opts.removeDups, "Remove duplicate Query-Target pairs (default: true)", true},
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
			RecordType:    fields[0],
			ClusterNumber: parseIntOrZero(fields[1]),
			Size:          parseIntOrZero(fields[2]),
			Identity:      parseFloatOrZero(fields[3]),
			Strand:        parseStrand(fields[4]),
			Unused1:       fields[5],
			Unused2:       fields[6],
			CIGAR:         fields[7],
			Query:         queryLabel,
			Target:        targetLabel,
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
		fmt.Fprintf(os.Stderr, "\033[31mucs: removed %d duplicate entries\033[0m\n", duplicateCount)
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
		_, err := fmt.Fprintf(writer, "recordType\tclusterNumber\tsize\tidentity\tstrand\tunused1\tunused2\tcigar\tquery\ttarget\n")
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
			if record.Identity != nil {
				identityStr = fmt.Sprintf("%.2f", *record.Identity)
			}

			_, err := fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
				record.RecordType,
				record.ClusterNumber,
				record.Size,
				identityStr,
				strandStr,
				record.Unused1,
				record.Unused2,
				record.CIGAR,
				record.Query,
				record.Target)
			if err != nil {
				return err
			}
		}
	}
	return nil
}

// Parquet writing function
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
			parquetRecords[i] = record.ToParquet()
		}

		// Write records
		if _, err := writer.Write(parquetRecords); err != nil {
			return fmt.Errorf("writing parquet records: %w", err)
		}
	}

	return nil
}

// Split sequence ID at semicolon
func splitSeqID(id string, split bool) string {
	if !split {
		return id
	}
	parts := strings.Split(id, ";")
	return parts[0]
}

// Parse integer fields (ClusterNumber and Size)
func parseIntOrZero(s string) uint32 {
	if s == "*" {
		return 0
	}
	val, err := strconv.ParseUint(s, 10, 32)
	if err != nil {
		return 0
	}
	return uint32(val)
}

// Parse float fields (Identity)
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

// Parse strand
func parseStrand(s string) *byte {
	if s == "+" {
		b := byte('+')
		return &b
	}
	if s == "-" {
		b := byte('-')
		return &b
	}
	return nil
}

// UC file summary
func summarizeUC(input *os.File, inputFileName string, opts Options) (int, int, int, int, error) {
	scanner, err := createScanner(input, inputFileName)
	if err != nil {
		return 0, 0, 0, 0, err
	}

	rowCount := 0
	querySequences := make(map[string]struct{})
	targetSequences := make(map[string]struct{})
	queryToTargets := make(map[string]map[string]struct{})

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

		// Initialize map for this query if it doesn't exist
		if _, exists := queryToTargets[queryLabel]; !exists {
			queryToTargets[queryLabel] = make(map[string]struct{})
		}

		// For centroid records (C), use Query as Target
		if fields[0] == "C" {
			targetSequences[queryLabel] = struct{}{}
			queryToTargets[queryLabel][queryLabel] = struct{}{}
		} else if fields[9] != "*" {
			targetLabel := splitSeqID(fields[9], opts.splitSeqID)
			targetSequences[targetLabel] = struct{}{}
			queryToTargets[queryLabel][targetLabel] = struct{}{}
		}
	}

	// Count queries mapped to multiple targets
	multiMappedQueries := 0
	for _, targets := range queryToTargets {
		if len(targets) > 1 {
			multiMappedQueries++
		}
	}

	if err := scanner.Err(); err != nil {
		return 0, 0, 0, 0, fmt.Errorf("reading input: %w", err)
	}

	return rowCount, len(querySequences), len(targetSequences), multiMappedQueries, nil
}

func writeSummary(output *os.File, rowCount, uniqueQuerySequences, uniqueTargetSequences, multiMappedQueries int) error {
	writer := bufio.NewWriter(output)
	defer writer.Flush()

	// Define the rows with their labels and values
	rows := []struct {
		label string
		value int
	}{
		{"Total lines in the file:", rowCount},
		{"Unique query sequences:", uniqueQuerySequences},
		{"Unique target sequences:", uniqueTargetSequences},
		{"Queries mapped to multiple targets:", multiMappedQueries},
	}

	// Find the longest label and the longest number
	maxLabelWidth := 0
	maxNumberWidth := 0
	for _, row := range rows {
		if len(row.label) > maxLabelWidth {
			maxLabelWidth = len(row.label)
		}
		numWidth := len(strconv.Itoa(row.value))
		if numWidth > maxNumberWidth {
			maxNumberWidth = numWidth
		}
	}

	// Create format string with calculated widths
	format := fmt.Sprintf("%%-%ds %%%dd\n", maxLabelWidth, maxNumberWidth)

	// Print each row
	for _, row := range rows {
		if row.label == "Queries mapped to multiple targets:" && row.value > 0 {
			_, err := fmt.Fprintf(writer, "\033[31m"+format+"\033[0m", row.label, row.value)
			if err != nil {
				return err
			}
		} else {
			_, err := fmt.Fprintf(writer, format, row.label, row.value)
			if err != nil {
				return err
			}
		}
	}

	return nil
}

// Buffered scanner for input file
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

// Convert UCRecord to ParquetRecord
func (r UCRecord) ToParquet() ParquetRecord {
	return ParquetRecord{
		RecordType:    r.RecordType,
		ClusterNumber: uint32(r.ClusterNumber),
		Size:          uint32(r.Size),
		Identity:      r.Identity,
		Strand:        r.Strand,
		Unused1:       r.Unused1,
		Unused2:       r.Unused2,
		CIGAR:         r.CIGAR,
		Query:         r.Query,
		Target:        r.Target,
	}
}
