# scSimEval Studio: Interactive Shiny Web Application

[![Posit Connect](https://img.shields.io/badge/Posit%20Connect-Cloud-blue?logo=rstudio)](https://kabilanbio-scsimeval.share.connect.posit.cloud/)
[![GitHub](https://img.shields.io/badge/GitHub-scSimEval-blue?logo=github)](https://github.com/kabilanbio/scSimEval)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**scSimEval Studio** is the dedicated, interactive Shiny web application for the [`scSimEval`](https://github.com/kabilanbio/scSimEval) R package. It provides an intuitive graphical interface for benchmarking single-cell transcriptomics (scRNA-seq), chromatin accessibility (scATAC-seq), and paired multiomics data simulation methods across 62 curated ground-truth-free evaluation measures.

---

## 🚀 Live Web Application

The application is hosted on **Posit Connect Cloud**:
👉 **[https://kabilanbio-scsimeval.share.connect.posit.cloud/](https://kabilanbio-scsimeval.share.connect.posit.cloud/)**

---

## 🛠️ Local Installation & Running

To run this Shiny application locally on your computer:

```r
# 1. Install required dependencies
install.packages(c("shiny", "bslib", "ggplot2", "DT", "Matrix", "remotes", "googleAnalyticsR"))

# 2. Install scSimEval package from GitHub
remotes::install_github("kabilanbio/scSimEval")

# 3. Launch the Shiny application
shiny::runApp()
```

---

## ☁️ Deployment to Posit Connect Cloud

This repository is pre-configured for automated deployment to Posit Connect Cloud:

```r
install.packages("rsconnect")
library(rsconnect)

rsconnect::deployApp(
  appDir  = ".",
  appName = "scsimeval",
  account = "kabilanbio"
)
```

### Environment Variables
Configure the following in the Posit Connect Cloud dashboard:
* `GA_MEASUREMENT_ID`: `G-D4BY0FVPTQ`
* `GA_PROPERTY_ID`: `557610038`
* `GA_AUTH_FILE`: `google_key.json`

---

## 👥 Authors & Maintainers

* **Sakthivel Kabilan** (Ph.D. Scholar, ICAR-IASRI, New Delhi) - [kabilan151414@gmail.com](mailto:kabilan151414@gmail.com)
* **Dr. Dwijesh Chandra Mishra** (Senior Scientist, ICAR-IASRI)
* **Dr. Shashi Bhushan Lal** (Principal Scientist, ICAR-IASRI)
* **Dr. Sudhir Srivastava** (Senior Scientist, ICAR-IASRI)
* **Dr. Krishna Kumar Chaturvedi** (Principal Scientist, ICAR-IASRI)
* **Dr. Sharanbasappa** (Scientist, ICAR-IASRI)

**Division of Agricultural Bioinformatics, ICAR-Indian Agricultural Statistics Research Institute (IASRI), New Delhi, India.**
