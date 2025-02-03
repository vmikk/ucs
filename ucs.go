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

// A type to store command options
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

// A type for Parquet output
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
		writer := bufio.NewWriter(output)
		defer writer.Flush()

		if strings.HasSuffix(opts.outputFile, ".parquet") {
			err = processAndWriteParquet(input, opts.outputFile, opts)
		} else {
			err = processAndWriteText(input, writer, opts)
		}
	}

	if err != nil {
		fatalError("Error processing file: %v", err)
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

// Streaming processor for text output
func processAndWriteText(input *os.File, writer *bufio.Writer, opts Options) error {
	scanner, err := createScanner(input, opts.inputFile)
	if err != nil {
		return err
	}

	seenPairs := make(map[string]struct{})
	queryToTargets := make(map[string]map[string]struct{})
	duplicateCount := 0

	// Write header first
	if opts.mapOnly {
		if _, err := writer.WriteString("Query\tTarget\n"); err != nil {
			return err
		}
	} else {
		if _, err := writer.WriteString("recordType\tclusterNumber\tsize\tidentity\tstrand\tunused1\tunused2\tcigar\tquery\ttarget\n"); err != nil {
			return err
		}
	}

	// Process records one by one
	for scanner.Scan() {
		record, ok := parseUCRecord(scanner.Text(), opts)
		if !ok {
			continue
		}

		if opts.removeDups {
			pairKey := record.Query + "\t" + record.Target
			if _, exists := seenPairs[pairKey]; exists {
				duplicateCount++
				continue
			}
			seenPairs[pairKey] = struct{}{}
		}

		if opts.multiMapped {
			if _, exists := queryToTargets[record.Query]; !exists {
				queryToTargets[record.Query] = make(map[string]struct{})
			}
			queryToTargets[record.Query][record.Target] = struct{}{}
		} else {
			if err := writeUCRecord(writer, record, opts); err != nil {
				return err
			}
		}
	}

	// Handle multi-mapped queries if needed
	if opts.multiMapped {
		for query, targets := range queryToTargets {
			if len(targets) > 1 {
				for target := range targets {
					record := UCRecord{Query: query, Target: target}
					if err := writeUCRecord(writer, record, opts); err != nil {
						return err
					}
				}
			}
		}
	}

	if duplicateCount > 0 {
		fmt.Fprintf(os.Stderr, "\033[31mucs: removed %d duplicate entries\033[0m\n", duplicateCount)
	}

	return scanner.Err()
}

// Streaming processor for Parquet output
func processAndWriteParquet(input *os.File, outputFile string, opts Options) error {
	f, err := os.Create(outputFile)
	if err != nil {
		return fmt.Errorf("creating output file: %w", err)
	}
	defer f.Close()

	scanner, err := createScanner(input, opts.inputFile)
	if err != nil {
		return err
	}

	// Configure ZSTD codec with better compression
	zstdCodec := &zstd.Codec{
		Level: zstd.SpeedBetterCompression,
	}

	if opts.mapOnly {
		type MapRecord struct {
			Query  string `parquet:"query"`
			Target string `parquet:"target"`
		}
		writer := parquet.NewGenericWriter[MapRecord](f, parquet.Compression(zstdCodec))
		defer writer.Close()

		seenPairs := make(map[string]struct{})
		for scanner.Scan() {
			record, ok := parseUCRecord(scanner.Text(), opts)
			if !ok {
				continue
			}

			if opts.removeDups {
				pairKey := record.Query + "\t" + record.Target
				if _, exists := seenPairs[pairKey]; exists {
					continue
				}
				seenPairs[pairKey] = struct{}{}
			}

			mapRecord := MapRecord{
				Query:  record.Query,
				Target: record.Target,
			}
			if _, err := writer.Write([]MapRecord{mapRecord}); err != nil {
				return fmt.Errorf("writing map record: %w", err)
			}
		}
	} else {
		writer := parquet.NewGenericWriter[ParquetRecord](f, parquet.Compression(zstdCodec))
		defer writer.Close()

		seenPairs := make(map[string]struct{})
		for scanner.Scan() {
			record, ok := parseUCRecord(scanner.Text(), opts)
			if !ok {
				continue
			}

			if opts.removeDups {
				pairKey := record.Query + "\t" + record.Target
				if _, exists := seenPairs[pairKey]; exists {
					continue
				}
				seenPairs[pairKey] = struct{}{}
			}

			parquetRecord := record.ToParquet()
			if _, err := writer.Write([]ParquetRecord{parquetRecord}); err != nil {
				return fmt.Errorf("writing parquet record: %w", err)
			}
		}
	}

	return scanner.Err()
}

// Helper function to write a single record
func writeUCRecord(writer *bufio.Writer, record UCRecord, opts Options) error {
	if opts.mapOnly {
		_, err := fmt.Fprintf(writer, "%s\t%s\n", record.Query, record.Target)
		return err
	}

	strandStr := "*"
	if record.Strand != nil {
		strandStr = string(*record.Strand)
	}

	identityStr := "*"
	if record.Identity != nil {
		identityStr = fmt.Sprintf("%.2f", *record.Identity)
	}

	_, err := fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
		record.RecordType, record.ClusterNumber, record.Size,
		identityStr, strandStr, record.Unused1, record.Unused2,
		record.CIGAR, record.Query, record.Target)
	return err
}

// UC file summary
func summarizeUC(input *os.File, inputFileName string, opts Options) (int, int, int, int, error) {
	scanner, err := createScanner(input, inputFileName)
	if err != nil {
		return 0, 0, 0, 0, err
	}

	rowCount := 0
	querySequences := make(map[string]struct{})            // Set of unique queries
	targetSequences := make(map[string]struct{})           // Set of unique targets
	queryToTargets := make(map[string]map[string]struct{}) // unique query to target pairs

	for scanner.Scan() {
		rowCount++ // Count every line
		line := scanner.Text()
		fields := strings.Split(line, "\t")

		// Skip broken lines
		if len(fields) < 10 {
			continue
		}

		// Skip C records for sequence analysis
		if fields[0] == "C" {
			continue
		}

		// Process all valid non-C records
		queryLabel := splitSeqID(fields[8], opts.splitSeqID)

		// Add query to the set of unique queries
		querySequences[queryLabel] = struct{}{}

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

// Parse a UC record from a line of text
func parseUCRecord(line string, opts Options) (UCRecord, bool) {
	fields := strings.Split(line, "\t")
	if len(fields) < 10 {
		return UCRecord{}, false
	}

	// Skip C records as they are redundant
	if fields[0] == "C" {
		return UCRecord{}, false
	}

	queryLabel := splitSeqID(fields[8], opts.splitSeqID)
	targetLabel := splitSeqID(fields[9], opts.splitSeqID)

	record := UCRecord{
		RecordType: fields[0],
		Query:      queryLabel,
		Target:     targetLabel,
	}

	// Process target based on record type
	switch record.RecordType {
	case "H":
		// Hit record - use target as is
		record.Target = targetLabel
	case "S":
		// Seed record - use query as both query and target
		record.Target = queryLabel
	case "N":
		// No hit - use query as target
		record.Target = queryLabel
	}

	// Parse optional fields for non-N and non-H records
	if record.RecordType != "H" && record.RecordType != "N" {
		if num, err := strconv.ParseUint(fields[1], 10, 32); err == nil {
			record.ClusterNumber = uint32(num)
		}
		if num, err := strconv.ParseUint(fields[2], 10, 32); err == nil {
			record.Size = uint32(num)
		}
		if fields[3] != "*" {
			if val, err := strconv.ParseFloat(fields[3], 64); err == nil {
				record.Identity = &val
			}
		}
		if fields[4] != "*" {
			strand := fields[4][0]
			record.Strand = &strand
		}
	}

	record.Unused1 = fields[5]
	record.Unused2 = fields[6]
	record.CIGAR = fields[7]

	return record, true
}

// Split sequence ID at semicolon if enabled
func splitSeqID(id string, split bool) string {
	if !split {
		return id
	}
	parts := strings.SplitN(id, ";", 2)
	return parts[0]
}
