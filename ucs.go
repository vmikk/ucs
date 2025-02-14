package main

import (
	"bufio"
	"compress/gzip"
	"flag"
	"fmt"
	"os"
	"strconv"
	"strings"
	"time"

	"github.com/briandowns/spinner"
	"github.com/parquet-go/parquet-go"
	"github.com/parquet-go/parquet-go/compress/zstd"
)

// Version information
const Version = "0.8.0"

// Custom error types for better error handling and testing
type UCError struct {
	Type    string
	Message string
	Err     error
}

func (e *UCError) Error() string {
	if e.Err != nil {
		return fmt.Sprintf("%s: %s: %v", e.Type, e.Message, e.Err)
	}
	return fmt.Sprintf("%s: %s", e.Type, e.Message)
}

func newUCError(errType, msg string, err error) *UCError {
	return &UCError{
		Type:    errType,
		Message: msg,
		Err:     err,
	}
}

// A type to store command options
type Options struct {
	inputFile   string
	outputFile  string
	summary     bool
	mapOnly     bool
	splitSeqID  bool
	removeDups  bool
	multiMapped bool
	version     bool
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

// A type for Parquet output (all columns)
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

// A type for simplified Parquet output
type MapRecord struct {
	Query  string `parquet:"query"`
	Target string `parquet:"target"`
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

// Helper function to check if we're in an interactive terminal session
func isInteractive() bool {
	// Check if stdin and stderr are terminals
	// We use stderr for spinner to avoid interfering with stdout redirections
	return isTerminal(os.Stdin) && isTerminal(os.Stderr)
}

// Helper function to check if file descriptor is a terminal
func isTerminal(f *os.File) bool {
	fileInfo, err := f.Stat()
	if err != nil {
		return false
	}
	return (fileInfo.Mode() & os.ModeCharDevice) != 0
}

const (
	green         = "\033[32m"
	reset         = "\033[0m"
	msgProcessing = green + "Processing UC data..." + reset
	msgComplete   = green + "Processing complete" + reset
)

// Create a spinner instance if we're in interactive mode
func createSpinner() *spinner.Spinner {
	if !isInteractive() {
		return nil
	}
	s := spinner.New(spinner.CharSets[14], 100*time.Millisecond,
		spinner.WithWriter(os.Stderr),
		spinner.WithSuffix(" "+msgProcessing))
	s.Color("green")
	return s
}

// Centralized error handling helper
func fatalError(format string, a ...interface{}) {
	// Red color code
	const red = "\033[31m"
	const reset = "\033[0m"

	msg := fmt.Sprintf(format, a...)
	fmt.Fprintf(os.Stderr, red+"Error: %s"+reset+"\n", msg)
	os.Exit(1)
}

func main() {
	opts := parseFlags()

	// Create and start spinner
	s := createSpinner()
	if s != nil {
		s.Start()
	}

	input, err := openInputFile(opts.inputFile)
	if err != nil {
		if s != nil {
			s.Stop()
		}
		fatalError("Error opening input file: %v", err)
	}
	defer input.Close()

	output, err := createOutputFile(opts.outputFile)
	if err != nil {
		if s != nil {
			s.Stop()
		}
		fatalError("Error creating output file: %v", err)
	}
	defer output.Close()

	if opts.summary {
		rowCount, uniqueQuerySequences, uniqueTargetSequences, multiMappedQueries, duplicateCount, err := summarizeUC(input, opts.inputFile, opts)
		if err != nil {
			if s != nil {
				s.Stop()
			}
			fatalError("Error processing file: %v", err)
		}
		// Stop spinner before writing summary
		if s != nil {
			s.Stop()
		}
		err = writeSummary(output, rowCount, uniqueQuerySequences, uniqueTargetSequences, multiMappedQueries, duplicateCount)
	} else {
		writer := bufio.NewWriter(output)
		defer writer.Flush()

		if strings.HasSuffix(opts.outputFile, ".parquet") {
			err = processAndWriteParquet(input, opts.outputFile, opts, s)
		} else {
			err = processAndWriteText(input, writer, opts, s)
		}
	}

	if err != nil {
		if s != nil {
			s.Stop()
		}
		fatalError("Error processing file: %v", err)
	}

	// Update spinner message on completion and stop
	if s != nil {
		s.FinalMSG = msgComplete + "\n"
		s.Stop()
	}
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
		{"version", "v", &opts.version, "Print version information", false},
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
		fmt.Fprintf(flag.CommandLine.Output(), `ucs %s - USEARCH/VSEARCH cluster format parser and converter

Usage:
  %s [flags]

Flags:
`, Version, os.Args[0])
		// Find the longest flag combination to determine padding
		maxLen := 0
		for _, f := range flagPairs {
			flagLen := len(f.long) + 2 // --flag
			if f.short != "" {
				flagLen += 4 // -x,
			}
			if flagLen > maxLen {
				maxLen = flagLen
			}
		}

		// Format string with consistent padding
		format := fmt.Sprintf("  %%-%ds\t%%s\n", maxLen)

		for _, f := range flagPairs {
			shortFlag := ""
			if f.short != "" {
				shortFlag = fmt.Sprintf("-%s, ", f.short)
			}
			flagName := fmt.Sprintf("%s--%s", shortFlag, f.long)
			fmt.Fprintf(flag.CommandLine.Output(), format, flagName, f.usage)
		}
		fmt.Fprintf(flag.CommandLine.Output(), "\nFor more information, visit https://github.com/vmikk/ucs\n")
	}

	flag.Parse()

	// Handle version flag early
	if opts.version {
		fmt.Printf("ucs %s\n", Version)
		os.Exit(0)
	}

	// Auto-detect stdin if no input file specified and stdin is a pipe
	if opts.inputFile == "-" && !isTerminal(os.Stdin) {
		opts.inputFile = "-"
	}

	// Check if any flags were provided or if stdin is a pipe
	if flag.NFlag() == 0 && len(flag.Args()) == 0 && isTerminal(os.Stdin) {
		fmt.Fprintf(os.Stderr, "\033[31mError: no arguments provided\033[0m\n\n")
		flag.Usage()
		os.Exit(1)
	}

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

// Split sequence ID at semicolon if enabled
func splitSeqID(id string, split bool) string {
	if !split {
		return id
	}
	parts := strings.SplitN(id, ";", 2)
	return parts[0]
}

// Parse the full UC record from a line of text
func parseUCRecord(line string, opts Options) (UCRecord, bool) {
	fields := strings.Split(line, "\t")

	// Skip broken lines
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
		// Parse cluster number and size for H records
		if num, err := strconv.ParseUint(fields[1], 10, 32); err == nil {
			record.ClusterNumber = uint32(num)
		}
		if num, err := strconv.ParseUint(fields[2], 10, 32); err == nil {
			record.Size = uint32(num)
		}
		// Parse identity and strand for H records
		if fields[3] != "*" {
			if val, err := strconv.ParseFloat(fields[3], 64); err == nil {
				record.Identity = &val
			}
		}
		if fields[4] != "*" {
			strand := fields[4][0]
			record.Strand = &strand
		}
	case "S":
		// Seed record - use query as both query and target
		record.Target = queryLabel
		// Parse cluster number and size for S records
		if num, err := strconv.ParseUint(fields[1], 10, 32); err == nil {
			record.ClusterNumber = uint32(num)
		}
		if num, err := strconv.ParseUint(fields[2], 10, 32); err == nil {
			record.Size = uint32(num)
		}
		// For S records, identity and strand are always "*"
	case "N":
		// No hit - use query as target
		record.Target = queryLabel
		// For N records, cluster and size should be parsed
		if num, err := strconv.ParseUint(fields[1], 10, 32); err == nil {
			record.ClusterNumber = uint32(num)
		}
		if num, err := strconv.ParseUint(fields[2], 10, 32); err == nil {
			record.Size = uint32(num)
		}
	}

	record.Unused1 = fields[5]
	record.Unused2 = fields[6]
	record.CIGAR = fields[7]

	return record, true
}

// Parse only Query and Target fields from a UC record
func parseMapRecord(line string, split bool) (UCRecord, bool) {
	// Split only up to field 10 (0-9)
	fields := strings.SplitN(line, "\t", 10)
	if len(fields) < 10 {
		return UCRecord{}, false
	}

	// Skip C records as they are redundant
	if fields[0] == "C" {
		return UCRecord{}, false
	}

	query := splitSeqID(fields[8], split)
	target := splitSeqID(fields[9], split)

	// Handle special cases based on record type
	switch fields[0] {
	case "S", "N":
		target = query
	}

	return UCRecord{
		RecordType: fields[0],
		Query:      query,
		Target:     target,
	}, true
}

// UC-file processing logic
func processRecords(scanner *bufio.Scanner, opts Options, handler func(UCRecord) error, s *spinner.Spinner) error {
	seenPairs := make(map[string]struct{})
	queryToTargets := make(map[string]map[string]struct{})
	duplicateCount := 0
	lineNum := 0

	for scanner.Scan() {
		lineNum++

		var err error
		if opts.mapOnly {
			// Optimized path for map-only mode
			record, ok := parseMapRecord(scanner.Text(), opts.splitSeqID)
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
				err = handler(record)
			}
		} else {
			// Original path for full record processing
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
				err = handler(record)
			}
		}

		if err != nil {
			return newUCError("IO", fmt.Sprintf("failed to write record at line %d", lineNum), err)
		}
	}

	// Handle multi-mapped queries if needed
	if opts.multiMapped {
		for query, targets := range queryToTargets {
			if len(targets) > 1 {
				for target := range targets {
					if err := handler(UCRecord{Query: query, Target: target}); err != nil {
						return newUCError("IO", fmt.Sprintf("failed to write multi-mapped record for query %s", query), err)
					}
				}
			}
		}
	}

	if duplicateCount > 0 {
		if s != nil {
			s.Stop() // Stop spinner before showing warning
			fmt.Fprintf(os.Stderr, "\033[31mucs: removed %d duplicate entries\033[0m\n", duplicateCount)
			s.Start() // Restart spinner for remaining processing
		} else {
			fmt.Fprintf(os.Stderr, "\033[31mucs: removed %d duplicate entries\033[0m\n", duplicateCount)
		}
	}

	return scanner.Err()
}

// Process UC-file and write output into TSV format
func processAndWriteText(input *os.File, writer *bufio.Writer, opts Options, s *spinner.Spinner) error {
	scanner, err := createScanner(input, opts.inputFile)
	if err != nil {
		return newUCError("IO", "failed to create scanner", err)
	}

	// Write header
	header := "Query\tTarget\n"
	if !opts.mapOnly {
		header = "recordType\tclusterNumber\tsize\tidentity\tstrand\tunused1\tunused2\tcigar\tquery\ttarget\n"
	}
	if _, err := writer.WriteString(header); err != nil {
		return newUCError("IO", "failed to write header", err)
	}

	return processRecords(scanner, opts, func(record UCRecord) error {
		return writeUCRecord(writer, record, opts)
	}, s)
}

// Process UC-file and write output into Parquet format
func processAndWriteParquet(input *os.File, outputFile string, opts Options, s *spinner.Spinner) error {
	f, err := os.Create(outputFile)
	if err != nil {
		return newUCError("IO", "failed to create output file", err)
	}
	defer f.Close()

	scanner, err := createScanner(input, opts.inputFile)
	if err != nil {
		return newUCError("IO", "failed to create scanner", err)
	}

	// Configure ZSTD codec with better compression
	zstdCodec := &zstd.Codec{Level: zstd.SpeedBetterCompression}

	if opts.mapOnly {
		writer := parquet.NewGenericWriter[MapRecord](f, parquet.Compression(zstdCodec))
		defer func() {
			if err := writer.Close(); err != nil {
				// Log the error since we can't return it from the defer
				fmt.Fprintf(os.Stderr, "\033[31mError closing parquet writer: %v\033[0m\n", err)
			}
		}()

		return processRecords(scanner, opts, func(record UCRecord) error {
			_, err := writer.Write([]MapRecord{{Query: record.Query, Target: record.Target}})
			return err
		}, s)
	}

	writer := parquet.NewGenericWriter[ParquetRecord](f, parquet.Compression(zstdCodec))
	defer func() {
		if err := writer.Close(); err != nil {
			fmt.Fprintf(os.Stderr, "\033[31mError closing parquet writer: %v\033[0m\n", err)
		}
	}()

	return processRecords(scanner, opts, func(record UCRecord) error {
		_, err := writer.Write([]ParquetRecord{record.ToParquet()})
		return err
	}, s)
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
func summarizeUC(input *os.File, inputFileName string, opts Options) (int, int, int, int, int, error) {
	scanner, err := createScanner(input, inputFileName)
	if err != nil {
		return 0, 0, 0, 0, 0, err
	}

	rowCount := 0
	querySequences := make(map[string]struct{})            // Set of unique queries
	targetSequences := make(map[string]struct{})           // Set of unique targets
	queryToTargets := make(map[string]map[string]struct{}) // Unique query to target pairs
	seenPairs := make(map[string]struct{})                 // Set to track duplicates
	duplicateCount := 0                                    // Number of duplicate query-target pairs

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
		targetLabel := splitSeqID(fields[9], opts.splitSeqID)
		if fields[0] == "S" {
			targetLabel = queryLabel
		}

		// Check for duplicates
		pairKey := queryLabel + "\t" + targetLabel
		if _, exists := seenPairs[pairKey]; exists {
			duplicateCount++
			continue
		}
		seenPairs[pairKey] = struct{}{}

		// Add query to the set of unique queries
		querySequences[queryLabel] = struct{}{}

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
		return 0, 0, 0, 0, 0, fmt.Errorf("reading input: %w", err)
	}

	return rowCount, len(querySequences), len(targetSequences), duplicateCount, multiMappedQueries, nil
}

func writeSummary(output *os.File, rowCount, uniqueQuerySequences, uniqueTargetSequences, duplicateCount, multiMappedQueries int) error {
	// Check if output is stdout
	useColors := output != os.Stdout

	// Define the rows with their labels and values
	rows := []struct {
		label string
		value int
		warn  bool // whether to print in red if value > 0
	}{
		{"Total lines in the file:", rowCount, false},
		{"Unique query sequences:", uniqueQuerySequences, false},
		{"Unique target sequences:", uniqueTargetSequences, false},
		{"Duplicate query-target pairs:", duplicateCount, true},
		{"Queries mapped to multiple targets:", multiMappedQueries, true},
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
	for _, row := range rows {
		if row.warn && row.value > 0 && useColors {
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
