package main

import (
	"testing"
	"time"

	. "github.com/onsi/ginkgo/v2"
	. "github.com/onsi/gomega"
)

func TestUcs(t *testing.T) {
	RegisterFailHandler(Fail)

	// Set global defaults for asynchronous assertions.
	SetDefaultEventuallyTimeout(10 * time.Second)
	SetDefaultEventuallyPollingInterval(100 * time.Millisecond)

	suiteConfig, reporterConfig := GinkgoConfiguration()
	// Use succinct mode for a clean, hierarchical output.
	reporterConfig.Verbose = false        // Disable verbose debug info.
	reporterConfig.Succinct = true        // Enable succinct output.
	reporterConfig.FullTrace = false      // Full stack traces not required.
	reporterConfig.ShowNodeEvents = false // Do not show BeforeEach/AfterEach events.
	reporterConfig.NoColor = false        // Enable colored output.

	RunSpecs(t, "UCS", suiteConfig, reporterConfig)
}

var _ = BeforeSuite(func() {
	By("Setting up test environment")
})

var _ = AfterSuite(func() {
	By("Cleaning up test environment")
})
