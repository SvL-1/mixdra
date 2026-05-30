# Code generator: ports the validated mixture-model functions from the Shiny
# app's model_functions.R into clean package files (R/models-binary.R and
# R/models-ternary.R). Run from the repo root:
#   R -q --vanilla -f tools/port_models.R
#
# The source functions already use explicit named arguments; porting is purely a
# rename of the argument names to the package's lower-case convention. To
# guarantee the maths is byte-identical we eval the source bodies verbatim and
# wrap them, then numerically cross-check every ported function against its
# source over random inputs (maxdiff must be 0).

src_path <- "MixTox_shiny_v2/functions/model_functions.R"
lines <- readLines(src_path, warn = FALSE)
# Strip the top-of-file library() calls (plotly etc. are not installed and are
# irrelevant to the pure-math functions); everything else is function defs.
lines <- lines[!grepl("^\\s*library\\(", lines)]
env <- new.env()
eval(parse(text = paste(lines, collapse = "\n")), envir = env)

# --- definitions of which source functions to port and how args map ----------
# new = package name; src = source name; bargs = base arg names in the package
# (in source order); extra = deviation arg names (package convention).
bin_base <- c("max", "slope1", "slope2", "ec501", "ec502")
tri_base <- c("max", "slope1", "slope2", "slope3", "ec50_1", "ec50_2", "ec50_3")

binary <- list(
  list(new="ca_bi",    src="CA_bi",    extra=character(0)),
  list(new="ia_bi",    src="IA_bi",    extra=character(0)),
  list(new="ca_sa_bi", src="CA_SA_bi", extra="a"),
  list(new="ca_dr_bi", src="CA_DR_bi", extra=c("a","b")),
  list(new="ca_dl_bi", src="CA_DL_bi", extra=c("a","b")),
  list(new="ia_sa_bi", src="IA_SA_bi", extra="a"),
  list(new="ia_dr_bi", src="IA_DR_bi", extra=c("a","b")),
  list(new="ia_dl_bi", src="IA_DL_bi", extra=c("a","b"))
)
ternary <- list(
  list(new="ca_tri",    src="CA",    extra=character(0)),
  list(new="ia_tri",    src="IA",    extra=character(0)),
  list(new="ca_sa_tri", src="CA_SA", extra="a"),
  list(new="ca_dr_tri", src="CA_DR", extra=c("a","b1","b2","b3")),
  list(new="ca_dl_tri", src="CA_DL", extra=c("a","b")),
  list(new="ia_sa_tri", src="IA_SA", extra="a"),
  list(new="ia_dr_tri", src="IA_DR", extra=c("a","b1","b2","b3")),
  list(new="ia_dl_tri", src="IA_DL", extra=c("a","b"))
)

# Verify our extra-arg assumptions against the actual source formals.
check_formals <- function(spec, base, nconc) {
  for (m in spec) {
    fa <- names(formals(get(m$src, envir = env)))
    expected_n <- nconc + length(base) + length(m$extra)
    if (length(fa) != expected_n)
      stop(sprintf("%s: source has %d args, expected %d (conc+base+extra)",
                   m$src, length(fa), expected_n))
  }
}
check_formals(binary, bin_base, 2)
check_formals(ternary, tri_base, 3)

# --- emit one package file ----------------------------------------------------
emit <- function(spec, base, conc, suffix_doc, path) {
  out <- c(
    sprintf("# %s mixture model functions.", suffix_doc),
    "#",
    "# Ported from MixTox_shiny_v2/functions/model_functions.R. Each source",
    "# function already takes explicit named arguments; the only change here is",
    "# renaming them to the package's lower-case convention. The bodies are",
    "# reproduced verbatim (via deparse of the source), so the maths is identical",
    "# (cross-checked numerically by tools/port_models.R: maxdiff == 0).",
    "")
  for (m in spec) {
    src_fn <- get(m$src, envir = env)
    new_args <- c(conc, base, m$extra)        # package arg names, source order
    old_args <- names(formals(src_fn))        # source arg names, same order
    body_txt <- deparse(body(src_fn))         # verbatim body, starts with "{"

    rox <- c(
      sprintf("#' %s: %s mixture predictor (ported verbatim from %s)",
              m$new, tolower(suffix_doc), m$src),
      "#'",
      sprintf("#' @param %s Concentrations.", paste(conc, collapse = ",")),
      sprintf("#' @param %s Model parameters.",
              paste(c(base, m$extra), collapse = ",")),
      "#' @return Predicted response (scalar).",
      "#' @keywords internal")

    header <- sprintf("%s <- function(%s) {", m$new, paste(new_args, collapse = ", "))
    # Re-pass the package args positionally into the verbatim source function,
    # which we re-create inline under its source signature.
    inner_hdr <- sprintf("  .fn <- function(%s) ", paste(old_args, collapse = ", "))
    inner <- c(paste0(inner_hdr, body_txt[1]), body_txt[-1])
    call_line <- sprintf("  .fn(%s)", paste(new_args, collapse = ", "))

    out <- c(out, rox, header, inner, call_line, "}", "")
  }
  # vectorised wrappers
  out <- c(out, sprintf("# Vectorised wrappers over the %s concentration vectors.",
                        paste(conc, collapse = "/")))
  va <- paste(sprintf('"%s"', conc), collapse = ", ")
  for (m in spec) {
    out <- c(out, sprintf("%s_vec <- Vectorize(%s, vectorize.args = c(%s))",
                          m$new, m$new, va))
  }
  out <- c(out, "")
  writeLines(out, path)
  cat("WROTE", path, "lines=", length(out), "\n")
}

emit(binary, bin_base, c("c1","c2"), "Binary (2-chemical)", "R/models-binary.R")
emit(ternary, tri_base, c("c1","c2","c3"), "Ternary (3-chemical)", "R/models-ternary.R")

# --- numerical cross-check: ported vs source over random inputs ---------------
# Load the freshly-written ports into their own env.
penv <- new.env()
eval(parse(text = paste(readLines("R/models-binary.R"), collapse="\n")), envir = penv)
eval(parse(text = paste(readLines("R/models-ternary.R"), collapse="\n")), envir = penv)

set.seed(123)
rand_conc <- function(n) sapply(seq_len(n), function(i) sample(c(0, runif(1, 0, 5)), 1))
rand_par  <- function(base, extra) {
  p <- c(runif(1, 100, 900))                              # max
  nslope <- sum(grepl("slope", base)); nec <- sum(grepl("ec50", base))
  p <- c(p, runif(nslope, 0.3, 6), runif(nec, 0.01, 3))   # slopes, ec50s
  if (length(extra)) p <- c(p, runif(length(extra), -5, 5))
  p
}
worst <- 0
report <- function(spec, base, nconc) {
  for (m in spec) {
    md <- 0
    for (i in 1:300) {
      cc <- rand_conc(nconc); pp <- rand_par(base, m$extra)
      a <- tryCatch(do.call(get(m$src, envir=env), c(as.list(cc), as.list(pp))),
                    error=function(e) NA)
      b <- tryCatch(do.call(get(m$new, envir=penv), c(as.list(cc), as.list(pp))),
                    error=function(e) NA)
      if (is.na(a) && is.na(b)) next
      if (xor(is.na(a), is.na(b))) { md <- Inf; next }
      md <- max(md, abs(a - b))
    }
    worst <<- max(worst, md)
    cat(sprintf("  %-10s vs %-8s maxdiff = %g\n", m$new, m$src, md))
  }
}
cat("Binary cross-check:\n");  report(binary, bin_base, 2)
cat("Ternary cross-check:\n"); report(ternary, tri_base, 3)
cat(if (worst == 0) "ALL_IDENTICAL\n" else sprintf("WORST=%g\n", worst))
