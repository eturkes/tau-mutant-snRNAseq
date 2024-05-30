#    This file is part of tau-mutant-snRNAseq.
#    Copyright (C) 2024  Emir Turkes, Naoto Watamura, UK DRI at UCL
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    Emir Turkes can be contacted at emir.turkes@eturkes.com

# This script runs all files in the analysis.
# Sections can be commented out as needed.

setwd(dirname(parent.frame(2)$ofile)) # Move to location of this file.

xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("misc", "allen_mouse_hip_ctx_10x.Rmd"),
    output_file = file.path(
      "..", "..", "results", "misc", "allen_mouse_hip_ctx_10x.html"
    ),
    envir = new.env()
  )
)

xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch01", "MAPTKI_batch01_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch01", "MAPTKI_batch01_01prep.html"
    ),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch01", "NLGF_MAPTKI_batch01_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch01", "NLGF_MAPTKI_batch01_01prep.html"
    ),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch01", "S305N_batch01_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch01", "S305N_batch01_01prep.html"
    ),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch01", "NLGF_S305N_batch01_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch01", "NLGF_S305N_batch01_01prep.html"
    ),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch01", "P301S_batch01_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch01", "P301S_batch01_01prep.html"
    ),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch01", "NLGF_P301S_batch01_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch01", "NLGF_P301S_batch01_01prep.html"
    ),
    envir = new.env()
  )
)

xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch02", "MAPTKI_batch02_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch02", "MAPTKI_batch02_01prep.html"
    ),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch02", "NLGF_MAPTKI_batch02_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch02", "NLGF_MAPTKI_batch02_01prep.html"
    ),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch02", "S305N_batch02_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch02", "S305N_batch02_01prep.html"
    ),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch02", "NLGF_S305N_batch02_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch02", "NLGF_S305N_batch02_01prep.html"
    ),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch02", "P301S_batch02_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch02", "P301S_batch02_01prep.html"
    ),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch02", "NLGF_P301S_batch02_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch02", "NLGF_P301S_batch02_01prep.html"
    ),
    envir = new.env()
  )
)

xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch03", "MAPTKI_batch03_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch03", "MAPTKI_batch03_01prep.html"
    ),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch03", "NLGF_MAPTKI_batch03_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch03", "NLGF_MAPTKI_batch03_01prep.html"
    ),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch03", "S305N_batch03_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch03", "S305N_batch03_01prep.html"
    ),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch03", "NLGF_S305N_batch03_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch03", "NLGF_S305N_batch03_01prep.html"
    ),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch03", "P301S_batch03_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch03", "P301S_batch03_01prep.html"
    ),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch03", "NLGF_P301S_batch03_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch03", "NLGF_P301S_batch03_01prep.html"
    ),
    envir = new.env()
  )
)

xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch04", "MAPTKI_batch04_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch04", "MAPTKI_batch04_01prep.html"
    ),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch04", "NLGF_MAPTKI_batch04_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch04", "NLGF_MAPTKI_batch04_01prep.html"
    ),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch04", "S305N_batch04_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch04", "S305N_batch04_01prep.html"
    ),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch04", "NLGF_S305N_batch04_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch04", "NLGF_S305N_batch04_01prep.html"
    ),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch04", "P301S_batch04_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch04", "P301S_batch04_01prep.html"
    ),
    envir = new.env()
  )
)
xfun::Rscript_call(
  rmarkdown::render,
  list(
    file.path("batch04", "NLGF_P301S_batch04_01prep.Rmd"),
    output_file = file.path(
      "..", "..", "results", "batch04", "NLGF_P301S_batch04_01prep.html"
    ),
    envir = new.env()
  )
)
