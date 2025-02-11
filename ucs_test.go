package main

import (
	// "bufio"
	// "compress/gzip"
	"os"
	// "path/filepath"
	//"strings"

	. "github.com/onsi/ginkgo/v2"
	. "github.com/onsi/gomega"
	//"github.com/parquet-go/parquet-go"
)

var _ = Describe("UCS", func() {
	const testFile = "test/test.uc.gz"
	var tmpDir string

	BeforeEach(func() {
		var err error
		tmpDir, err = os.MkdirTemp("", "ucs-test-*")
		Expect(err).NotTo(HaveOccurred())
	})

	AfterEach(func() {
		os.RemoveAll(tmpDir)
	})

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

})
