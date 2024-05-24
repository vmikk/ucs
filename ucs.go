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
	inputFile := flag.String("input", "", "Input file (default: stdin)")
	outputFile := flag.String("output", "", "Output file (default: stdout)")
	flag.StringVar(inputFile, "i", "", "Input file (default: stdin)")
	flag.StringVar(outputFile, "o", "", "Output file (default: stdout)")
	flag.Usage = func() {
		fmt.Fprintf(flag.CommandLine.Output(), "Usage of %s:\n", os.Args[0])
		fmt.Fprintf(flag.CommandLine.Output(), "  -i, --input  <file/stdin>\tInput  file (default: stdin)\n")
		fmt.Fprintf(flag.CommandLine.Output(), "  -o, --output <file/stdout>\tOutput file (default: stdout)\n")
		// flag.PrintDefaults()
	}
	flag.Parse()

	var input *os.File
	var output *os.File
	var err error

	if *inputFile == "-" {
		input = os.Stdin
	} else {
		input, err = os.Open(*inputFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
			os.Exit(1)
		}
		defer input.Close()
	}

	if *outputFile == "-" {
		output = os.Stdout
	} else {
		output, err = os.Create(*outputFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
			os.Exit(1)
		}
		defer output.Close()
	}

	var scanner *bufio.Scanner

	// Check if the input file is gzip-compressed
	if strings.HasSuffix(*inputFile, ".gz") {
		gzipReader, err := gzip.NewReader(input)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating gzip reader: %v\n", err)
			os.Exit(1)
		}
		defer gzipReader.Close()
		scanner = bufio.NewScanner(gzipReader)
	} else {
		scanner = bufio.NewScanner(input)
	}

	rowCount := 0
	querySequences := make(map[string]struct{})
	targetSequences := make(map[string]struct{})

	// Process records
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
		fmt.Fprintf(os.Stderr, "Error reading input: %v\n", err)
		os.Exit(1)
	}

	writer := bufio.NewWriter(output)
	fmt.Fprintf(writer, "RowCount\t%d\n", rowCount)
	fmt.Fprintf(writer, "UniqueQuerySequences\t%d\n", len(querySequences))
	fmt.Fprintf(writer, "UniqueTargetSequences\t%d\n", len(targetSequences))
	writer.Flush()
}
