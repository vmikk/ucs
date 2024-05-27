package main

import (
	"bufio"
	"compress/gzip"
	"flag"
	"fmt"
	"os"
	"strings"
)

func main() {
	inputFile, outputFile := parseFlags()

	input, err := openInputFile(*inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
		os.Exit(1)
	}
	defer input.Close()

	output, err := createOutputFile(*outputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
		os.Exit(1)
	}
	defer output.Close()

	rowCount, uniqueQuerySequences, uniqueTargetSequences, err := processFile(input, *inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error processing file: %v\n", err)
		os.Exit(1)
	}

	err = writeOutput(output, rowCount, uniqueQuerySequences, uniqueTargetSequences)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error writing output: %v\n", err)
		os.Exit(1)
	}
}

func parseFlags() (*string, *string) {
	inputFile := flag.String("input", "-", "Input file (default: stdin)")
	outputFile := flag.String("output", "-", "Output file (default: stdout)")
	flag.StringVar(inputFile, "i", "-", "Input file (default: stdin)")
	flag.StringVar(outputFile, "o", "-", "Output file (default: stdout)")
	flag.Usage = func() {
		fmt.Fprintf(flag.CommandLine.Output(), "Usage of %s:\n", os.Args[0])
		fmt.Fprintf(flag.CommandLine.Output(), "  -i, --input  <file/stdin>\tInput   file (default: stdin)\n")
		fmt.Fprintf(flag.CommandLine.Output(), "  -o, --output <file/stdout>\tOutput file (default: stdout)\n")
	}
	flag.Parse()
	return inputFile, outputFile
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

func processFile(input *os.File, inputFileName string) (int, int, int, error) {
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
		rowCount++
		line := scanner.Text()
		fields := strings.Split(line, "\t")
		if len(fields) < 10 {
			continue
		}

		query := fields[8]
		target := fields[9]

		querySequences[query] = struct{}{}
		if target != "*" {
			targetSequences[target] = struct{}{}
		} else {
			targetSequences[query] = struct{}{}
		}
	}

	if err := scanner.Err(); err != nil {
		return 0, 0, 0, fmt.Errorf("reading input: %w", err)
	}

	return rowCount, len(querySequences), len(targetSequences), nil
}

func writeOutput(output *os.File, rowCount, uniqueQuerySequences, uniqueTargetSequences int) error {
	writer := bufio.NewWriter(output)
	defer writer.Flush()
	_, err := fmt.Fprintf(writer, "RowCount\t%d\n", rowCount)
	if err != nil {
		return err
	}
	_, err = fmt.Fprintf(writer, "UniqueQuerySequences\t%d\n", uniqueQuerySequences)
	if err != nil {
		return err
	}
	_, err = fmt.Fprintf(writer, "UniqueTargetSequences\t%d\n", uniqueTargetSequences)
	return err
}
