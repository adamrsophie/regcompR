# Export tables and plots

`export_table()` writes a
[`regcomp_table()`](https://adamrsophie.github.io/regcompR/reference/regcomp_table.md)
(or the table of one or several fits) to a file whose format follows the
extension: `.csv`, `.html`, `.tex` (LaTeX), `.md` (Markdown) or `.docx`
(Word, with the flextable package). `export_plot()` draws a diagnostic
plot of a fit (see
[`plot.regcomp_fit()`](https://adamrsophie.github.io/regcompR/reference/plot.regcomp_fit.md))
into a `.png`, `.jpeg`, `.tiff`, `.pdf` or `.svg` file.

## Usage

``` r
export_table(x, file, ...)

export_plot(fit, file, which = "qq", width = 7, height = 5, res = 300, ...)
```

## Arguments

- x:

  A `regcomp_table`, a fitted model, or a list of fitted models.

- file:

  Output file.

- ...:

  For `export_table()`, arguments to
  [`regcomp_table()`](https://adamrsophie.github.io/regcompR/reference/regcomp_table.md)
  when `x` holds fits; for `export_plot()`, arguments to
  [`plot.regcomp_fit()`](https://adamrsophie.github.io/regcompR/reference/plot.regcomp_fit.md).

- fit:

  A fitted model.

- which:

  The plot to draw.

- width, height:

  Size in inches.

- res:

  Resolution in pixels per inch for bitmap formats.

## Value

The file name, invisibly.

## Examples

``` r
d <- data.frame(x = rnorm(100)); d$y <- rpois(100, exp(0.5 + 0.3 * d$x))
f <- regcomp(y ~ x, family = "poisson", data = d)
export_table(f, tempfile(fileext = ".csv"))
export_plot(f, tempfile(fileext = ".png"), which = "qq")
```
