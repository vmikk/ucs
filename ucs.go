package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
)

func main() {
	inputFile := flag.String("input", "", "Input file (default: stdin)")
	outputFile := flag.String("output", "", "Output file (default: stdout)")
	flag.Parse()

	var input *os.File
	var output *os.File
	var err error

	if *inputFile == "" {
		input = os.Stdin
	} else {
		input, err = os.Open(*inputFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
			os.Exit(1)
		}
		defer input.Close()
	}

	if *outputFile == "" {
		output = os.Stdout
	} else {
		output, err = os.Create(*outputFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
			os.Exit(1)
		}
		defer output.Close()
	}

	scanner := bufio.NewScanner(input)
	rowCount := 0

	for scanner.Scan() {
		rowCount++
	}

	if err := scanner.Err(); err != nil {
		fmt.Fprintf(os.Stderr, "Error reading input: %v\n", err)
		os.Exit(1)
	}

	writer := bufio.NewWriter(output)
	fmt.Fprintf(writer, "RowCount\t%d\n", rowCount)
	writer.Flush()
}
