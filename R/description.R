# Step 1: Update DESCRIPTION file
writeLines(c(
  "Package: R2DAG",
  "Type: Package",
  "Title: Used codes for the R2 connectedness article",
  "Version: 0.1.0",
  "Authors@R: person(\"Marton\", \"Espan\", email = \"espan.marci@email.com\", role = c(\"aut\", \"cre\"))",
  "Description: A package for R2 connectedness analysis.",
  "License: MIT",
  "Encoding: UTF-8",
  "LazyData: true",
  "RoxygenNote: 7.1.1",
  "Imports:",
  "    pcalg,",
  "    igraph,",
  "    relaimpo"
), "DESCRIPTION")

# Step 2: Update NAMESPACE file (Regenerate if using roxygen2)
writeLines(c(
  "#' @import pcalg",
  "#' @import igraph",
  "#' @import relaimpo",
  "NULL"
), "NAMESPACE")

# Step 3: Document the package (regenerates NAMESPACE)
devtools::document()
