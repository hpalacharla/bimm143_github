# Class 18: Pertussis and the CMI-PB Project
Harshitha Palacharla (PID: A17775151)

- [Section 1. Investigating pertussis cases by
  year](#section-1-investigating-pertussis-cases-by-year)
- [Section 2. A tale of two vaccines (wP &
  aP)](#section-2-a-tale-of-two-vaccines-wp--ap)
- [Section 3. Exploring CMI-PB data](#section-3-exploring-cmi-pb-data)
  - [Side-Note: Working with dates](#side-note-working-with-dates)
  - [Joining multiple tables](#joining-multiple-tables)
- [Section 4. Examine IgG Ab titer
  levels](#section-4-examine-igg-ab-titer-levels)
- [Section 5. Obtaining CMI-PB RNASeq
  data](#section-5-obtaining-cmi-pb-rnaseq-data)

## Section 1. Investigating pertussis cases by year

> Q1. With the help of the R “addin” package datapasta assign the CDC
> pertussis case number data to a data frame called `cdc` and use
> **ggplot** to make a plot of cases numbers over time.

Note: Using `read.csv` instead, upon instruction that it is okay to do
so in place of datapasta.

``` r
# Read in cdc pertussis data
cdc <- read.csv("U.S. Reported Pertussis Cases_ 1922 - 2025.csv")
```

``` r
# Convert pertussis case numbers from strings to numerics

cdc$Number.of.Reported.Pertussis.Cases <- as.numeric(gsub("," , "" ,
                                                          cdc$Number.of.Reported.Pertussis.Cases))
```

``` r
# Plot case numbers over time

library(ggplot2)

ggplot(cdc) +
  aes(x = Year, y = Number.of.Reported.Pertussis.Cases) +
  geom_point() +
  geom_line() +
  labs(title = "Pertussis Cases by Year (1922-2019)", x = "Year", y = "Number of Cases")
```

![](class18_files/figure-commonmark/unnamed-chunk-3-1.png)

## Section 2. A tale of two vaccines (wP & aP)

> Q2. Using the ggplot geom_vline() function add lines to your previous
> plot for the 1946 introduction of the wP vaccine and the 1996 switch
> to aP vaccine (see example in the hint below). What do you notice?

``` r
# Add lines for wP introduction and switch to aP

ggplot(cdc) +
  aes(x = Year, y = Number.of.Reported.Pertussis.Cases) +
  geom_point() +
  geom_line() +
  labs(title = "Pertussis Cases by Year (1922-2019)", x = "Year", y = "Number of Cases") +
  geom_vline(xintercept = 1946, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 1996, color = "red", linetype = "dashed") +
  geom_text(x = 1946, y = 250000, label = "wP", color = "blue", hjust = -0.1) +
  geom_text(x = 1996, y = 250000, label = "aP", color = "red", hjust = -0.1)
```

![](class18_files/figure-commonmark/unnamed-chunk-4-1.png)

It appears that there was a substantial decline in pertussis cases
following the introduction of the wP vaccine, in 1946. Though levels of
pertussis cases initially remained steadily low past 1975, a shift
occurred following the introduction of the aP vaccine: cases are shown
to rise again following 1996.

> Q3. Describe what happened after the introduction of the aP vaccine?
> Do you have a possible explanation for the observed trend?

Pertussis cases seem to have been *increasing* since the introduction of
the aP vaccine in 1996. This is in contrast to the downward trend of
cases after the original, wP vaccine’s introduction.

It is possible that the increase in cases is due to the aP vaccine
granting a shorter period of protecion/immunity from pertussis, as
compared to the wP vaccine. Thus, the immunity of younger generations
vaccinated with aP may be reducing/waning sooner following childhood
vaccination, resulting in increased susceptibility to pertussis and more
pertussis cases. That being said, the following are also possible:
diagnosis methods have improved (i.e. better or more accessible
testing), fewer people are getting vaccinated (i.e. due to vaccine
hesistancy, reduced caution, etc.), or the pertussis bacteria has
evolved.

## Section 3. Exploring CMI-PB data

``` r
# Load jsonlite package to read, write and process JSON data
library(jsonlite)
```

``` r
# Read in metadata about individuals in study
subject <- read_json("https://www.cmi-pb.org/api/subject", simplifyVector = TRUE) 

# Take a quick look at the metadata
head(subject, 3)
```

      subject_id infancy_vac biological_sex              ethnicity  race
    1          1          wP         Female Not Hispanic or Latino White
    2          2          wP         Female Not Hispanic or Latino White
    3          3          wP         Female                Unknown White
      year_of_birth date_of_boost      dataset
    1    1986-01-01    2016-09-12 2020_dataset
    2    1968-01-01    2019-01-28 2020_dataset
    3    1983-01-01    2016-10-10 2020_dataset

> Q4. How many aP and wP infancy vaccinated subjects are in the dataset?

``` r
table(subject$infancy_vac)
```


    aP wP 
    87 85 

There are 87 aP vaccinated and 85 wP vaccinated subjects.

> Q5. How many Male and Female subjects/patients are in the dataset?

``` r
table(subject$biological_sex)
```


    Female   Male 
       112     60 

There are 60 Male subjects and 112 Female subjects.

> Q6. What is the breakdown of race and biological sex (e.g. number of
> Asian females, White males etc…)?

``` r
table(subject$biological_sex, subject$race)
```

            
             American Indian/Alaska Native Asian Black or African American
      Female                             0    32                         2
      Male                               1    12                         3
            
             More Than One Race Native Hawaiian or Other Pacific Islander
      Female                 15                                         1
      Male                    4                                         1
            
             Unknown or Not Reported White
      Female                      14    48
      Male                         7    32

The breakdown by race and biological sex is shown in the table above. It
appears the most represented groups are White females (48), followed by
White males and Asian females (these last two groups are tied in count,
32).

There is comparatively much less representation of the American
Indian/Alaska Native, Black/African American, and Pacific Islander
races. There is also higher representation of females than males.

### Side-Note: Working with dates

``` r
# Example

library(lubridate)
```

``` r
# Today's date
today()
```

    [1] "2026-06-08"

``` r
# Days since new year 2000
today() - ymd("2000-01-01")
```

    Time difference of 9655 days

``` r
# Convert above answer to years
time_length(today() - ymd("2000-01-01"), "years")
```

    [1] 26.43395

Today is “2026-05-28”. 9644 days have passed since new year 2000, which
is 26.4 years.

> Q7. Using this approach determine (i) the average age of wP
> individuals, (ii) the average age of aP individuals; and (iii) are
> they significantly different?

``` r
# Store age (in days) of subjects
subject$age <- today() - ymd(subject$year_of_birth)
```

``` r
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
# Filter to wP-vaccinated subjects.
wp <- subject |>
  filter(infancy_vac == "wP")

# Examine 6 number summary of their ages.
round ( summary (time_length(wp$age, "years")) )
```

       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
         23      33      35      37      40      58 

1)  The average age of wP individuals is 37 years.

``` r
# Filter to aP-vaccinated subjects.
ap <- subject |>
  filter(infancy_vac == "aP")

# Examine 6 number summary of their ages.
round ( summary (time_length(ap$age, "years")) )
```

       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
         23      27      28      28      29      35 

2)  The average age of aP individuals is 28 years.

``` r
# Run t.test of wP and aP ages to determine if there's a stat difference
t.test(time_length(wp$age, "years"), time_length(ap$age, "years"))
```


        Welch Two Sample t-test

    data:  time_length(wp$age, "years") and time_length(ap$age, "years")
    t = 12.918, df = 104.03, p-value < 2.2e-16
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
      7.407351 10.094058
    sample estimates:
    mean of x mean of y 
      37.0800   28.3293 

``` r
ggplot(subject) +
  aes(x = as.factor(infancy_vac), y = time_length(age, "years")) +
geom_boxplot() +
  xlab("Vaccination type") +
  ylab("Age in years")
```

![](class18_files/figure-commonmark/unnamed-chunk-16-1.png)

3)  Results from a two-sample t test of the ages of wP-vaccinated
    subjects and aP-vaccinated subjects indicate that the difference in
    ages *is statistically significant*. The p-value is very low
    (2.2e-16), so we find evidence for the alternative hypothesis that
    the true difference in mean ages is not equal to 0. On average, the
    wP individuals (mean age = 37 years) are older than the aP
    individuals (mean age = 28 years).

Even from the summary stats and a representative boxplot, we observe
that the middle 50% of the wP group (i.e. 1st quartile age to 3rd
quartile age) is 33-40 years old, whereas the middle 50% of the aP group
is younger at 27-29 years old.

> Q8. Determine the age of all individuals at time of boost?

``` r
age_at_boost <- time_length( ymd(subject$date_of_boost) - ymd(subject$year_of_birth), "years")

head(age_at_boost)
```

    [1] 30.69678 51.07461 33.77413 28.65982 25.65914 28.77481

> Q9. With the help of a faceted boxplot or histogram (see below), do
> you think these two groups are significantly different?

``` r
ggplot(subject) +
  aes(time_length(age, "year"),
      fill=as.factor(infancy_vac)) +
  geom_histogram(show.legend=FALSE) +
  facet_wrap(vars(infancy_vac), nrow=2) +
  xlab("Age in years")
```

    `stat_bin()` using `bins = 30`. Pick better value `binwidth`.

![](class18_files/figure-commonmark/unnamed-chunk-18-1.png)

The two groups appear to have a significant difference in age: the
wP-vaccinated group is older than the aP-vaccinated group. The
distribution of ages for the wP group is shifted right on the histogram
(toward higher ages), relative to the distribution for the aP group.
Furthermore, the t-test conducted earlier found evidence for a
difference, with a p-value of 2.2e-16.

### Joining multiple tables

``` r
# Read the specimen and ab_titer tables into R
specimen <- read_json("https://www.cmi-pb.org/api/specimen", simplifyVector = TRUE) 
titer <- read_json("https://www.cmi-pb.org/api/plasma_ab_titer", simplifyVector = TRUE) 
```

> Q9. Complete the code to join `specimen` and `subject` tables to make
> a new merged data frame containing all specimen records along with
> their associated subject details:

``` r
# Merge specimen and subject tables.
meta <- inner_join(specimen, subject)
```

    Joining with `by = join_by(subject_id)`

``` r
dim(meta)
```

    [1] 1503   14

``` r
head(meta)
```

      specimen_id subject_id actual_day_relative_to_boost
    1           1          1                           -3
    2           2          1                            1
    3           3          1                            3
    4           4          1                            7
    5           5          1                           11
    6           6          1                           32
      planned_day_relative_to_boost specimen_type visit infancy_vac biological_sex
    1                             0         Blood     1          wP         Female
    2                             1         Blood     2          wP         Female
    3                             3         Blood     3          wP         Female
    4                             7         Blood     4          wP         Female
    5                            14         Blood     5          wP         Female
    6                            30         Blood     6          wP         Female
                   ethnicity  race year_of_birth date_of_boost      dataset
    1 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    2 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    3 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    4 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    5 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    6 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
             age
    1 14768 days
    2 14768 days
    3 14768 days
    4 14768 days
    5 14768 days
    6 14768 days

- Note: `inner_join` was used, to ensure all rows have both a
  specimen_id and an associated subject_id.

> Q10. Now using the same procedure join `meta` with `titer` data so we
> can further analyze this data in terms of time of visit aP/wP,
> male/female etc.

``` r
abdata <- inner_join(titer, meta)
```

    Joining with `by = join_by(specimen_id)`

``` r
dim(abdata)
```

    [1] 52576    21

> Q11. How many specimens (i.e. entries in abdata) do we have for each
> isotype?

``` r
table(abdata$isotype)
```


      IgE   IgG  IgG1  IgG2  IgG3  IgG4 
     6698  5389 10117 10124 10124 10124 

There are 6698 specimens of IgE, 5389 specimens of IgG, 10117 specimens
of IgG1, and 10124 specimens each of IgG2, IgG3, and IgG4.

> Q12. What are the different \$dataset values in abdata and what do you
> notice about the number of rows for the most “recent” dataset?

``` r
table(abdata$dataset)
```


    2020_dataset 2021_dataset 2022_dataset 2023_dataset 
           31520         8085         7301         5670 

`$dataset` has 4 values, representing the 4 independent datasets/cohorts
created over 4 years: 2020_dataset, 2021_dataset, 2022_dataset,
2023_dataset. By creating a new, independent dataset each year, the
researchers could use the prior years’ datasets as “training data” for
their computational prediction models, and subsequently test the models
on the newest dataset.

The most recent dataset, from 2023, has the fewest number of rows (5670,
compared to 31520 in the first, 2020 dataset).

## Section 4. Examine IgG Ab titer levels

``` r
# Filter abdata for IgG isotype
igg <- abdata %>% filter(isotype == "IgG")
head(igg)
```

      specimen_id isotype is_antigen_specific antigen        MFI MFI_normalised
    1           1     IgG                TRUE      PT   68.56614       3.736992
    2           1     IgG                TRUE     PRN  332.12718       2.602350
    3           1     IgG                TRUE     FHA 1887.12263      34.050956
    4          19     IgG                TRUE      PT   20.11607       1.096366
    5          19     IgG                TRUE     PRN  976.67419       7.652635
    6          19     IgG                TRUE     FHA   60.76626       1.096457
       unit lower_limit_of_detection subject_id actual_day_relative_to_boost
    1 IU/ML                 0.530000          1                           -3
    2 IU/ML                 6.205949          1                           -3
    3 IU/ML                 4.679535          1                           -3
    4 IU/ML                 0.530000          3                           -3
    5 IU/ML                 6.205949          3                           -3
    6 IU/ML                 4.679535          3                           -3
      planned_day_relative_to_boost specimen_type visit infancy_vac biological_sex
    1                             0         Blood     1          wP         Female
    2                             0         Blood     1          wP         Female
    3                             0         Blood     1          wP         Female
    4                             0         Blood     1          wP         Female
    5                             0         Blood     1          wP         Female
    6                             0         Blood     1          wP         Female
                   ethnicity  race year_of_birth date_of_boost      dataset
    1 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    2 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    3 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    4                Unknown White    1983-01-01    2016-10-10 2020_dataset
    5                Unknown White    1983-01-01    2016-10-10 2020_dataset
    6                Unknown White    1983-01-01    2016-10-10 2020_dataset
             age
    1 14768 days
    2 14768 days
    3 14768 days
    4 15864 days
    5 15864 days
    6 15864 days

> Q13. Complete the following code to make a summary boxplot of Ab titer
> levels (MFI) for all antigens.

``` r
ggplot(igg) +
  aes(MFI_normalised, antigen) +
  geom_boxplot() + 
    xlim(0,75) +
  facet_wrap(vars(visit), nrow=2)
```

    Warning: Removed 5 rows containing non-finite outside the scale range
    (`stat_boxplot()`).

![](class18_files/figure-commonmark/unnamed-chunk-25-1.png)

> Q14. What antigens show differences in the level of IgG antibody
> titers recognizing them over time? Why these and not others?

The antigens associated with differences over time include PT (pertussis
toxin), PRN (pertactin autotransporter),FIM2/3 (Fimbrial protein 2/3),
and FHA (Filamentous hemagglutinin). This makes sense since the booster
vaccines use purified acellular pertussis antigens from the *Bordetella
pertussis* bacterium–including PT, PRN, FIM2/3, and FHA–to trigger an
immune response and pertussis-targeted antibody production.

The OVA (ovalbumin) associated IgG antibody level does not increase
since OVA is not present in the vaccine formulation, instead being used
as a negative control in the assay. DT (diphtheria toxoid) and TT
(tetanus toxin), though also present in the booster vaccines along with
the aP antigens, do not seem to show significant change in their
associated IgG antibody titer levels. This may be due, in part to the
reduced dosage of DT in the Tdap booster; additionally, the initial
childhood vaccination may offer better long-term immunity/recognition
against DT/TT than against pertussis, resulting in the Tdap booster
triggering a stronger immune response to the pertussis-specific
antigens.

``` r
# Make plot showing differences between the wP and aP groups

ggplot(igg) +
  aes(MFI_normalised, antigen, col=infancy_vac ) +
  geom_boxplot(show.legend = FALSE) + 
  facet_wrap(vars(visit), nrow=2) +
  xlim(0,75) +
  theme_bw()
```

    Warning: Removed 5 rows containing non-finite outside the scale range
    (`stat_boxplot()`).

![](class18_files/figure-commonmark/unnamed-chunk-26-1.png)

``` r
# Add infancy-vac to the faceting, for an easier-to-read plot

igg %>% filter(visit != 8) %>%
ggplot() +
  aes(MFI_normalised, antigen, col=infancy_vac ) +
  geom_boxplot(show.legend = FALSE) + 
  xlim(0,75) +
  facet_wrap(vars(infancy_vac, visit), nrow=2)
```

    Warning: Removed 5 rows containing non-finite outside the scale range
    (`stat_boxplot()`).

![](class18_files/figure-commonmark/unnamed-chunk-27-1.png)

> Q15. Filter to pull out only two specific antigens for analysis and
> create a boxplot for each. You can chose any you like. And the same
> for FIM2/3.

Boxplot for OVA (negative control):

``` r
filter(igg, antigen=="OVA") %>%
  ggplot() +
  aes(MFI_normalised, col=infancy_vac) +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(vars(visit)) +
  theme_bw() +
  labs(title = "OVA Antigen Levels per Visit (wP group in teal, aP group in red)")
```

![](class18_files/figure-commonmark/unnamed-chunk-28-1.png)

Boxplot for PT (pertussis toxin, an aP antigen):

``` r
filter(igg, antigen=="PT") %>%
  ggplot() +
  aes(MFI_normalised, col=infancy_vac) +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(vars(visit)) +
  theme_bw() +
  labs(title = "PT Antigen Levels per Visit (wP group in teal, aP group in red)")
```

![](class18_files/figure-commonmark/unnamed-chunk-29-1.png)

Boxplot for FIM2/3 (Fimbrial protein 2/3, an aP antigen)

``` r
filter(igg, antigen=="FIM2/3") %>%
  ggplot() +
  aes(MFI_normalised, col=infancy_vac) +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(vars(visit)) +
  theme_bw() +
  labs(title = "FIM2/3 Antigen Levels per Visit (wP group in teal, aP group in red)")
```

![](class18_files/figure-commonmark/unnamed-chunk-30-1.png)

Note: The following, *optional* plots were generated for OVA, PT, and
FIM2/3, respectively:

``` r
filter(igg, antigen=="OVA") %>%
  ggplot() +
    aes(as.factor(visit), MFI_normalised) +
    geom_boxplot(show.legend = FALSE) +
    facet_wrap(vars(infancy_vac)) +
    theme_bw() +
  labs(title = "OVA antigen levels over time, post-Tdap injection")
```

![](class18_files/figure-commonmark/unnamed-chunk-31-1.png)

``` r
filter(igg, antigen=="PT") %>%
  ggplot() +
    aes(as.factor(visit), MFI_normalised) +
    geom_boxplot(show.legend = FALSE) +
    facet_wrap(vars(infancy_vac)) +
    theme_bw() +
  labs(title = "PT antigen levels over time, post-Tdap injection")
```

![](class18_files/figure-commonmark/unnamed-chunk-32-1.png)

``` r
filter(igg, antigen=="FIM2/3") %>%
  ggplot() +
    aes(as.factor(visit), MFI_normalised) +
    geom_boxplot(show.legend = FALSE) +
    facet_wrap(vars(infancy_vac)) +
    theme_bw() +
  labs(title = "FIM2/3 antigen levels over time, post-Tdap injection")
```

![](class18_files/figure-commonmark/unnamed-chunk-33-1.png)

> Q16. What do you notice about these two antigens time courses and the
> PT data in particular?

The OVA (negative control) levels do not increase over time for either
the aP or wP groups, staying relatively constant post-Tdap injection. In
contrast, the PT (pertussis toxin; an aP antigen) levels are elevated by
visit 4, reach a peak at visit 5, and diminish significantly by visit 8.
Following a similar trend, the FIM2/3 (another aP antigen) levels are
also elevated by visit 4, though they peak at visit 7 and stay elevated
(i.e. do not diminish significantly afterwards).

> Q17. Do you see any clear difference in aP vs. wP responses?

There does not appear to be any clear difference in the aP versus wP
responses, both of which show very similar changes over time (as
described in Q16) in PT and FIM2/3 IgG titer levels. Likewise, neither
group shows an increase in OVA levels over time, as expected for the
negative control.

Note: It is worth mentioning that the PT IgG titer levels appear to
reach higher levels in any given visit for the wP group than the aP
group. However, this is specific to the PT antigen and is not observed
for FIM2/3.

``` r
# Plot time course data for IgG PT antigen levels, from the 2021 dataset

abdata.21 <- abdata %>% filter(dataset == "2021_dataset")

abdata.21 %>% 
  filter(isotype == "IgG",  antigen == "PT") %>%
  ggplot() +
    aes(x=planned_day_relative_to_boost,
        y=MFI_normalised,
        col=infancy_vac,
        group=subject_id) +
    geom_point() +
    geom_line() +
    geom_vline(xintercept=0, linetype="dashed") +
    geom_vline(xintercept=14, linetype="dashed") +
  labs(title="2021 dataset IgG PT",
       subtitle = "Dashed lines indicate days 0 (pre-boost) and 14 (apparent peak levels)")
```

![](class18_files/figure-commonmark/unnamed-chunk-34-1.png)

> Q18. Does this trend look similar for the 2020 dataset?

``` r
# Plot time course data for IgG PT antigen levels, from the 2020 dataset

abdata.20 <- abdata %>% filter(dataset == "2020_dataset")

abdata.20 %>% 
  filter(isotype == "IgG",  antigen == "PT") %>%
  ggplot() +
    aes(x=planned_day_relative_to_boost,
        y=MFI_normalised,
        col=infancy_vac,
        group=subject_id) +
    geom_point() +
    geom_line() +
    geom_vline(xintercept=0, linetype="dashed") +
    geom_vline(xintercept=14, linetype="dashed") +
  labs(title="2020 dataset IgG PT",
       subtitle = "Dashed lines indicate days 0 (pre-boost) and 14")
```

![](class18_files/figure-commonmark/unnamed-chunk-35-1.png)

``` r
# The same plot as above for 2020, taking a closer look at the first 125 days relative to boost

abdata.20 <- abdata %>% filter(dataset == "2020_dataset")

abdata.20 %>% 
  filter(isotype == "IgG",  antigen == "PT") %>%
  ggplot() +
    aes(x=planned_day_relative_to_boost,
        y=MFI_normalised,
        col=infancy_vac,
        group=subject_id) +
    geom_point() +
    geom_line() +
  coord_cartesian(xlim = c(0, 125)) +
    geom_vline(xintercept=0, linetype="dashed") +
    geom_vline(xintercept=14, linetype="dashed") +
  labs(title="2020 dataset IgG PT, zoomed in",
       subtitle = "Dashed lines indicate days 0 (pre-boost) and 14")
```

![](class18_files/figure-commonmark/unnamed-chunk-36-1.png)

The 2020 dataset does appear to show the same overall trend as the 2021
dataset. In both datasets, IgG PT titer levels seem to generally peak at
day 14 relative to receiving the booster injection, and decrease in the
days following.

It is worth noting, however, that the 2020 dataset shows more
variability than the 2021 dataset (perhaps due to having more data, as
seen in Q12), with some peaking–particularly among the wP-vaccinated
group–occuring at day 30 post-boost instead. However, the peak IgG PT
titer level still, more generally, remains at day 14.

## Section 5. Obtaining CMI-PB RNASeq data

``` r
# Read in available RNA-Seq data for the IGHG1 gene
url <- "https://www.cmi-pb.org/api/v2/rnaseq?versioned_ensembl_gene_id=eq.ENSG00000211896.7"

rna <- read_json(url, simplifyVector = TRUE) 
```

``` r
# Join the rna expression data with our metadata meta (i.e. sample + specimen data).
ssrna <- inner_join(rna, meta)
```

    Joining with `by = join_by(specimen_id)`

> Q19. Make a plot of the time course of gene expression for IGHG1 gene
> (i.e. a plot of visit vs. tpm).

``` r
ggplot(ssrna) +
  aes(visit, tpm, group=subject_id) +
  geom_point() +
  geom_line(alpha=0.2)
```

![](class18_files/figure-commonmark/unnamed-chunk-39-1.png)

> Q20.: What do you notice about the expression of this gene (i.e. when
> is it at its maximum level)?

The expression of the IGHG1 gene generally appears to peak at visit 4.
There is variability, with some peaking occuring at later visits;
however, the overall trend remains that maximum gene expression largely
occurs at visit 4.

> Q21. Does this pattern in time match the trend of antibody titer data?
> If not, why not?

Yes, the pattern in time does match. The trend of antibody titer data
for IgG PT, taken as an example, peaked around visit 5 (see Q15, 16). It
makes sense, then, that peak IGHG1 gene expression generally occurs
beforehand, at visit 4: since IGHG1 expression drives production of the
IgG1 antibody, it is reasonable and expected for maximum gene expression
to precede maximum antibody titer levels.

``` r
# Color and facet IGHG1 gene expression time course by infancy_vac status

ggplot(ssrna) +
  aes(tpm, col=infancy_vac) +
  geom_boxplot() +
  facet_wrap(vars(visit))
```

![](class18_files/figure-commonmark/unnamed-chunk-40-1.png)

``` r
# Density plot of IGHG1 gene expression, focusing on visit 4 (no clear wP vs. aP difference)

ssrna %>%  
  filter(visit==4) %>% 
  ggplot() +
    aes(tpm, col=infancy_vac) + geom_density() + 
    geom_rug() 
```

![](class18_files/figure-commonmark/unnamed-chunk-41-1.png)
