package main

import (
	"bufio"
	"compress/gzip"
	"flag"
	"fmt"
	"os"
	"strconv"
	"strings"
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
	ClusterNumber          int
	SeqLengthOrClusterSize int
	PercentIdentity        float64
	Strand                 string
	UnusedField1           string
	UnusedField2           string
	Alignment              string
	QueryLabel             string
	TargetLabel            string
	Query                  string // Processed query name
	OTU                    string // Processed OTU name
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
		rowCount, uniqueQuerySequences, uniqueTargetSequences, err := summarizeUC(input, opts.inputFile)
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

	flag.StringVar(&opts.inputFile, "input", "-", "Input file (default: stdin)")
	flag.StringVar(&opts.outputFile, "output", "-", "Output file (default: stdout)")
	flag.StringVar(&opts.inputFile, "i", "-", "Input file (default: stdin)")
	flag.StringVar(&opts.outputFile, "o", "-", "Output file (default: stdout)")
	flag.BoolVar(&opts.summary, "summary", false, "Print summary statistics (default: false)")
	flag.BoolVar(&opts.summary, "s", false, "Print summary statistics (default: false)")
	flag.BoolVar(&opts.mapOnly, "map-only", true, "Output only Query-OTU mapping (default: true)")
	flag.BoolVar(&opts.mapOnly, "m", true, "Output only Query-OTU mapping (default: true)")
	flag.BoolVar(&opts.splitSeqID, "split-id", true, "Split sequence IDs at semicolon (default: true)")
	flag.BoolVar(&opts.removeDups, "rm-dups", true, "Remove duplicate entries (default: true)")

	flag.Usage = func() {
		fmt.Fprintf(flag.CommandLine.Output(), "Usage of %s:\n", os.Args[0])
		fmt.Fprintf(flag.CommandLine.Output(), "  -i, --input    <file/stdin>\tInput file (default: stdin)\n")
		fmt.Fprintf(flag.CommandLine.Output(), "  -o, --output   <file/stdout>\tOutput file (default: stdout)\n")
		fmt.Fprintf(flag.CommandLine.Output(), "  -s, --summary\t\t\tPrint summary statistics (default: false)\n")
		fmt.Fprintf(flag.CommandLine.Output(), "  -m, --map-only\t\tOutput only Query-OTU mapping (default: true)\n")
		fmt.Fprintf(flag.CommandLine.Output(), "      --split-id\t\tSplit sequence IDs at semicolon (default: true)\n")
		fmt.Fprintf(flag.CommandLine.Output(), "      --rm-dups\t\t\tRemove duplicate entries (default: true)\n")
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
		// Create new scanner that includes the peeked bytes
		scanner = bufio.NewScanner(reader)
	}

	if isGzipped {
		gzipReader, err := gzip.NewReader(input)
		if err != nil {
			return nil, fmt.Errorf("creating gzip reader: %w", err)
		}
		defer gzipReader.Close()
		scanner = bufio.NewScanner(gzipReader)
	} else if scanner == nil { // if scanner wasn't created for stdin
		scanner = bufio.NewScanner(input)
	}

	var records []UCRecord
	seenPairs := make(map[string]struct{}) // For duplicate detection

	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.Split(line, "\t")
		if len(fields) < 10 {
			continue
		}

		// Skip S records - they are redundant with C records
		if fields[0] == "S" {
			continue
		}

		// Parse record
		record := UCRecord{
			RecordType:             fields[0],
			ClusterNumber:          parseIntOrZero(fields[1]),
			SeqLengthOrClusterSize: parseIntOrZero(fields[2]),
			PercentIdentity:        parseFloatOrZero(fields[3]),
			Strand:                 fields[4],
			UnusedField1:           fields[5],
			UnusedField2:           fields[6],
			Alignment:              fields[7],
			QueryLabel:             fields[8],
			TargetLabel:            fields[9],
		}

		// Process Query name
		record.Query = splitSeqID(record.QueryLabel, opts.splitSeqID)

		// Process OTU name based on record type
		switch record.RecordType {
		case "H":
			// Hit record - OTU comes from target label
			record.OTU = splitSeqID(record.TargetLabel, opts.splitSeqID)
		case "C":
			// Cluster centroid - OTU is the query sequence
			record.OTU = record.Query
		default:
			// Skip other record types (N records)
			continue
		}

		// Check for duplicates if removal is enabled
		if opts.removeDups {
			pairKey := record.Query + "\t" + record.OTU
			if _, exists := seenPairs[pairKey]; exists {
				continue
			}
			seenPairs[pairKey] = struct{}{}
		}

		records = append(records, record)
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("reading input: %w", err)
	}

	return records, nil
}

func writeUCRecords(writer *bufio.Writer, records []UCRecord, opts Options) error {
	if opts.mapOnly {
		// Write header for Query-OTU mapping
		_, err := fmt.Fprintf(writer, "Query\tOTU\n")
		if err != nil {
			return err
		}

		// Write data
		for _, record := range records {
			_, err := fmt.Fprintf(writer, "%s\t%s\n", record.Query, record.OTU)
			if err != nil {
				return err
			}
		}
	} else {
		// Write header for full table
		_, err := fmt.Fprintf(writer, "recordType\tclusterNumber\tseqLengthOrClusterSize\tpercentIdentity\tstrand\tunusedField1\tunusedField2\talignment\tqueryLabel\ttargetLabel\n")
		if err != nil {
			return err
		}

		// Write data
		for _, record := range records {
			_, err := fmt.Fprintf(writer, "%s\t%d\t%d\t%g\t%s\t%s\t%s\t%s\t%s\t%s\n",
				record.RecordType,
				record.ClusterNumber,
				record.SeqLengthOrClusterSize,
				record.PercentIdentity,
				record.Strand,
				record.UnusedField1,
				record.UnusedField2,
				record.Alignment,
				record.Query,
				record.OTU)
			if err != nil {
				return err
			}
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

func parseIntOrZero(s string) int {
	val, err := strconv.Atoi(s)
	if err != nil {
		return 0
	}
	return val
}

func parseFloatOrZero(s string) float64 {
	if s == "*" {
		return 0
	}
	val, err := strconv.ParseFloat(s, 64)
	if err != nil {
		return 0
	}
	return val
}

func summarizeUC(input *os.File, inputFileName string) (int, int, int, error) {
	var scanner *bufio.Scanner
	if strings.HasSuffix(inputFileName, ".gz") {
		gzipReader, err := gzip.NewReader(input)
		if err != nil {
			return 0, 0, 0, fmt.Errorf("creating gzip reader: %w", err)
		}
		defer gzipReader.Close()
		scanner = bufio.NewScanner(gzipReader)
	} else {
		scanner = bufio.NewScanner(input)
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
		querySequences[fields[8]] = struct{}{}

		// For centroid records (C), use Query as Target
		if fields[0] == "C" {
			targetSequences[fields[8]] = struct{}{}
		} else if fields[9] != "*" {
			targetSequences[fields[9]] = struct{}{}
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
