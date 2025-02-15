package main

import (
	"bufio"
	// "compress/gzip"
	"os"
	"path/filepath"
	"strings"

	. "github.com/onsi/ginkgo/v2"
	. "github.com/onsi/gomega"
	"github.com/parquet-go/parquet-go"
)

var _ = Describe("UCS", func() {

	// Test file
	const testFile = "test/test.uc.gz"

	var tmpDir string

	// Create a temporary directory for each test
	BeforeEach(func() {
		var err error
		tmpDir, err = os.MkdirTemp("", "ucs-test-*")
		Expect(err).NotTo(HaveOccurred())
	})

	// Clean up the temporary directory after each test
	AfterEach(func() {
		os.RemoveAll(tmpDir)
	})

	// ---------- Summary mode ----------

	Context("Summary mode", func() {
		It("should generate correct statistics for test.uc.gz", func() {
			input, err := openInputFile("./test/test.uc.gz")
			Expect(err).NotTo(HaveOccurred())
			defer input.Close()

			opts := Options{
				splitSeqID:  true,
				removeDups:  true,
				multiMapped: false,
			}

			rowCount, uniqueQueries, uniqueTargets, dups, multiMapped, err := summarizeUC(input, "test/test.uc.gz", opts)
			Expect(err).NotTo(HaveOccurred())

			// Assert known values
			Expect(rowCount).To(Equal(25329))
			Expect(uniqueQueries).To(Equal(24953))
			Expect(uniqueTargets).To(Equal(376))
			Expect(dups).To(Equal(0))
			Expect(multiMapped).To(Equal(0))
		})
	})

	// ---------- Map-only mode ----------

	Context("Map-only mode", func() {
		It("should correctly process TSV output with default settings", func() {
			outFile := filepath.Join(tmpDir, "out.tsv")

			// Run ucs with map-only mode
			opts := Options{
				inputFile:   testFile,
				outputFile:  outFile,
				mapOnly:     true,
				splitSeqID:  true,
				removeDups:  true,
				multiMapped: false,
			}

			// Process the file
			input, err := openInputFile(opts.inputFile)
			Expect(err).NotTo(HaveOccurred())
			defer input.Close()

			output, err := createOutputFile(opts.outputFile)
			Expect(err).NotTo(HaveOccurred())
			defer output.Close()

			writer := bufio.NewWriter(output)
			err = processAndWriteText(input, writer, opts, nil)
			Expect(err).NotTo(HaveOccurred())
			writer.Flush()

			// Read and verify output
			content, err := os.ReadFile(outFile)
			Expect(err).NotTo(HaveOccurred())

			lines := strings.Split(string(content), "\n")

			// Check header
			Expect(lines[0]).To(Equal("Query\tTarget"))

			// Check number of columns
			for _, line := range lines[1:] {
				if line == "" {
					continue
				}
				cols := strings.Split(line, "\t")
				Expect(cols).To(HaveLen(2))
			}

			// Count unique queries and targets
			queries := make(map[string]struct{})
			targets := make(map[string]struct{})
			for _, line := range lines[1:] {
				if line == "" {
					continue
				}
				cols := strings.Split(line, "\t")
				queries[cols[0]] = struct{}{}
				targets[cols[1]] = struct{}{}
			}

			// Verify expected counts
			Expect(len(queries)).To(BeNumerically("==", 24953))
			Expect(len(targets)).To(BeNumerically("==", 376))
		})

		It("should correctly process Parquet output", func() {
			outFile := filepath.Join(tmpDir, "out.parquet")

			opts := Options{
				inputFile:   testFile,
				outputFile:  outFile,
				mapOnly:     true,
				splitSeqID:  true,
				removeDups:  true,
				multiMapped: false,
			}

			// Process the file
			input, err := openInputFile(opts.inputFile)
			Expect(err).NotTo(HaveOccurred())
			defer input.Close()

			err = processAndWriteParquet(input, outFile, opts, nil)
			Expect(err).NotTo(HaveOccurred())

			// Read and verify Parquet output
			f, err := os.Open(outFile)
			Expect(err).NotTo(HaveOccurred())
			defer f.Close()

			reader := parquet.NewGenericReader[MapRecord](f)
			defer reader.Close()

			records := make([]MapRecord, reader.NumRows())
			n, err := reader.Read(records)
			Expect(err).NotTo(HaveOccurred())
			Expect(n).To(Equal(len(records)))

			// Verify records
			queries := make(map[string]struct{})
			targets := make(map[string]struct{})
			for _, record := range records {
				queries[record.Query] = struct{}{}
				targets[record.Target] = struct{}{}
			}

			Expect(len(queries)).To(BeNumerically("==", 24953))
			Expect(len(targets)).To(BeNumerically("==", 376))
		})
	})

})
