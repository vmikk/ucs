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
	inputFile   string
	outputFile  string
	summary     bool
	mapOnly     bool
	splitSeqID  bool
	removeDups  bool
	multiMapped bool
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
	Strand        string   `parquet:"strand"`
	Unused1       string   `parquet:"unused_1"`
	Unused2       string   `parquet:"unused_2"`
	CIGAR         string   `parquet:"cigar"`
	Query         string   `parquet:"query"`
	Target        string   `parquet:"target"`
}

// Convert UCRecord to ParquetRecord
func (r UCRecord) ToParquet() ParquetRecord {
	// Convert strand byte pointer to string
	var strandStr string
	if r.Strand != nil {
		strandStr = string(*r.Strand)
	} else {
		strandStr = "*"
	}

	return ParquetRecord{
		RecordType:    r.RecordType,
		ClusterNumber: r.ClusterNumber,
		Size:          r.Size,
		Identity:      r.Identity,
		Strand:        strandStr,
		Unused1:       r.Unused1,
		Unused2:       r.Unused2,
		CIGAR:         r.CIGAR,
		Query:         r.Query,
		Target:        r.Target,
	}
}

func main() {
	opts := parseFlags()

	input, err := openInputFile(opts.inputFile)
	if err != nil {
		fatalError("Error opening input file: %v", err)
	}
	defer input.Close()

	output, err := createOutputFile(opts.outputFile)
	if err != nil {
		fatalError("Error creating output file: %v", err)
	}
	defer output.Close()

	if opts.summary {
		rowCount, uniqueQuerySequences, uniqueTargetSequences, multiMappedQueries, err := summarizeUC(input, opts.inputFile, opts)
		if err != nil {
			fatalError("Error processing file: %v", err)
		}
		err = writeSummary(output, rowCount, uniqueQuerySequences, uniqueTargetSequences, multiMappedQueries)
	} else {
		records, err := processUCFile(input, opts.inputFile, opts)
		if err != nil {
			fatalError("Error processing file: %v", err)
		}

		writer := bufio.NewWriter(output)
		err = writeUCRecords(writer, records, opts)
		if err != nil {
			fatalError("Error writing output: %v", err)
		}
		writer.Flush()
	}

	if err != nil {
		fatalError("Error writing output: %v", err)
	}
}

// Centralized error handling helper
func fatalError(format string, a ...interface{}) {
	fmt.Fprintf(os.Stderr, format+"\n", a...)
	os.Exit(1)
}

// Parse command line flags
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
		{"split-id", "S", &opts.splitSeqID, "Split sequence IDs at semicolon (default: true)", true},
		{"rm-dups", "d", &opts.removeDups, "Remove duplicate Query-Target pairs (default: true)", true},
		{"multi-mapped", "M", &opts.multiMapped, "Output only queries mapped to multiple targets", false},
	}

	// Register all flags using a helper function
	registerFlag := func(long, short string, value interface{}, usage string, def interface{}) {
		switch v := value.(type) {
		case *string:
			flag.StringVar(v, long, def.(string), usage)
			if short != "" {
				flag.StringVar(v, short, def.(string), usage)
			}
		case *bool:
			flag.BoolVar(v, long, def.(bool), usage)
			if short != "" {
				flag.BoolVar(v, short, def.(bool), usage)
			}
		}
	}

	for _, f := range flagPairs {
		registerFlag(f.long, f.short, f.value, f.usage, f.def)
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
	seenPairs := make(map[string]struct{})        // For duplicate detection
	queryToTargets := make(map[string][]UCRecord) // Track targets per query
	duplicateCount := 0

	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.Split(line, "\t")

		// Skip broken lines
		if len(fields) < 10 {
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
		case "S":
			// Seed record - use query as both query and target
			record.Target = queryLabel
		case "C", "N":
			// Skip C records (redundant with S records) and N records (no hits)
			continue
		}

		// Check for duplicates if removal is enabled
		if opts.removeDups {
			var pairKeyBuilder strings.Builder
			pairKeyBuilder.WriteString(record.Query)
			pairKeyBuilder.WriteByte('\t')
			pairKeyBuilder.WriteString(record.Target)
			pairKey := pairKeyBuilder.String()

			if _, exists := seenPairs[pairKey]; exists {
				duplicateCount++
				continue
			}
			seenPairs[pairKey] = struct{}{}
		}

		if opts.multiMapped {
			queryToTargets[record.Query] = append(queryToTargets[record.Query], record)
		} else {
			records = append(records, record)
		}
	}

	// If multi-mapped mode is enabled, filter for queries with multiple targets
	if opts.multiMapped {
		for _, queryRecords := range queryToTargets {
			// Get unique targets for this query
			targets := make(map[string]struct{})
			for _, r := range queryRecords {
				targets[r.Target] = struct{}{}
			}

			// If query maps to multiple targets, add all its records
			if len(targets) > 1 {
				records = append(records, queryRecords...)
			}
		}
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

	var sb strings.Builder
	if opts.mapOnly {
		// Write header for Query-OTU mapping
		sb.WriteString("Query\tOTU\n")

		// Write data
		for _, record := range records {
			sb.WriteString(record.Query)
			sb.WriteByte('\t')
			sb.WriteString(record.Target)
			sb.WriteByte('\n')
		}
	} else {
		// Write header for full table
		sb.WriteString("recordType\tclusterNumber\tsize\tidentity\tstrand\tunused1\tunused2\tcigar\tquery\ttarget\n")

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

			sb.WriteString(record.RecordType)
			sb.WriteByte('\t')
			sb.WriteString(strconv.FormatUint(uint64(record.ClusterNumber), 10))
			sb.WriteByte('\t')
			sb.WriteString(strconv.FormatUint(uint64(record.Size), 10))
			sb.WriteByte('\t')
			sb.WriteString(identityStr)
			sb.WriteByte('\t')
			sb.WriteString(strandStr)
			sb.WriteByte('\t')
			sb.WriteString(record.Unused1)
			sb.WriteByte('\t')
			sb.WriteString(record.Unused2)
			sb.WriteByte('\t')
			sb.WriteString(record.CIGAR)
			sb.WriteByte('\t')
			sb.WriteString(record.Query)
			sb.WriteByte('\t')
			sb.WriteString(record.Target)
			sb.WriteByte('\n')
		}
	}

	_, err := writer.WriteString(sb.String())
	return err
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

		// Skip broken lines
		// Skip C records as they are redundant
		if len(fields) < 10 || fields[0] == "C" {
			continue
		}

		rowCount++
		queryLabel := splitSeqID(fields[8], opts.splitSeqID)

		targetLabel := splitSeqID(fields[9], opts.splitSeqID)
		if fields[0] == "S" {
			targetLabel = queryLabel
		}

		if _, exists := queryToTargets[queryLabel]; !exists {
			queryToTargets[queryLabel] = make(map[string]struct{})
		}

		if targetLabel != "*" {
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

	var sb strings.Builder
	// Print each row
	for _, row := range rows {
		if row.label == "Queries mapped to multiple targets:" && row.value > 0 {
			sb.WriteString(fmt.Sprintf("\033[31m"+format+"\033[0m", row.label, row.value))
		} else {
			sb.WriteString(fmt.Sprintf(format, row.label, row.value))
		}
	}

	writer := bufio.NewWriter(output)
	defer writer.Flush()
	_, err := writer.WriteString(sb.String())
	return err
}

// Buffered scanner for input file
func createScanner(input *os.File, inputFileName string) (*bufio.Scanner, error) {
	reader := bufio.NewReader(input)

	// Check if input is gzipped, either by filename or content
	isGzipped := strings.HasSuffix(inputFileName, ".gz")
	if !isGzipped && inputFileName == "-" {
		// For stdin, peek at the first two bytes to check for gzip magic number
		magic, err := reader.Peek(2)
		if err == nil && len(magic) == 2 && magic[0] == 0x1f && magic[1] == 0x8b {
			isGzipped = true
		}
	}

	if isGzipped {
		gzipReader, err := gzip.NewReader(reader)
		if err != nil {
			return nil, fmt.Errorf("creating gzip reader: %w", err)
		}
		return bufio.NewScanner(gzipReader), nil
	}

	return bufio.NewScanner(reader), nil
}
