#!/bin/bash

set -e  # Exit immediately if a command fails

# Step 1: Install R dependencies
Rscript "./requirements.R"

# Step 2: Compile Rcpp attributes
echo "Running Rcpp::compileAttributes()..."
Rscript -e "Rcpp::compileAttributes()"

# Step 3: Generate documentation
echo "Running roxygen2::roxygenise()..."
Rscript -e "roxygen2::roxygenise()"

# Step 4: Get package name and version from DESCRIPTION
PKG_NAME=$(grep -E "^Package:" DESCRIPTION | awk '{print $2}')
PKG_VERSION=$(grep -E "^Version:" DESCRIPTION | awk '{print $2}')
TARBALL="${PKG_NAME}_${PKG_VERSION}.tar.gz"

# Step 5: Build tarball
echo "Building tarball..."
R CMD build .

# Step 6: Run checks
echo "Running R CMD check..."
R CMD check "$TARBALL"

# Step 7: Install the package
echo "Installing the package..."
R CMD INSTALL --preclean .

echo "Done!"
