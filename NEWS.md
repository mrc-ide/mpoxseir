# mpoxseir 0.2.11
* Can fit to cumulative cases by age with nested binomials

# mpoxseir 0.2.10

* Now allow "burundi" as a option for region
* Bug fixed where there was no seeding in ASW in Sud-Kivu

# mpoxseir 0.2.9

* children_ind_raw replaced with is_child, adults_ind_raw replaced with (1 - is_child)
* vaccination allocation determined based on 3 criteria of eligibility, prioritisation and target met

# mpoxseir 0.2.8

* Age groups are now 0-4, 5-11, 12-14, 15-19, then in 5yr bands to allow for new guidance on vaccines for +/-12yo
* Add optional synthetic contact matrices from Prem et al. 2021

# mpoxseir 0.2.7

* `beta_hcw` scaled by total population size

# mpoxseir 0.2.6

* Those aged 12 and above are considered adults (relevant to 12 - 17 & CSWs)

# mpoxseir 0.2.5

* Fixed lag issue affecting some states

# mpoxseir 0.2.4

* Increase population of SWs in Sud Kivu

# mpoxseir 0.2.3

* Export cases by transmission route

# mpoxseir 0.2.2

* Added fitting to cumulative CFR by age group

# mpoxseir 0.2.1

* Export cases and deaths in key populations and vaccine doses by age group / key pop

# mpoxseir 0.2.0

* Package ported from odin-dust-mcstate to odin2-dust2-monty

# mpoxseir 0.1.11

* Parametrised HCWs
* Added fitting to % cases in HCW / SW
* Fixed output bug in cases among SW

# mpoxseir 0.1.10

* Added Negative-Binomial likelihood option

# mpoxseir 0.1.7

* Redefined mpoxseir date so that "2023-01-01" corresponds to mpoxseir date 0

# mpoxseir 0.1.6

* Add mpoxseir_date functions

# mpoxseir 0.1.4

* Split SW into CSW and ASW, added HCW

# mpoxseir 0.0.1

* Basic age-structured SEIR in odin.dust, adapted from squire
