package main_test

import (
	"testing"

	. "github.com/onsi/ginkgo/v2"
	. "github.com/onsi/gomega"
)

func TestUcs(t *testing.T) {
	RegisterFailHandler(Fail)
	RunSpecs(t, "ucs test suite")
}
