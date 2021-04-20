# Research Log
*Alain Chauffoureaux*

#### 20 April 2021
##### Removing supplementary information from Web of Life networks:
In many networks from Web of Life, there are columns/rows that are used to give supplementary information, usually species abundance.
For now, I simply remove rows/columns that match predefined names and do not use this additional information.
There are most likely many supplementary rows/columns that I did not remove yet.
Since these rows will not be matched against taxonomic databases, I will remove them manually when I correct species names.

##### Cleaning and resolving species names:
Today I start cleaning and resolving species names.
I have read the supplementary methods from F. Cagua that explain in details how he resolved species names in his study.
Her is a short list of the steps he mentions:
- [ ] Listing species from networks
- [ ] Determine hierarchical level (subspecies, species, genus) and aggregate at species level
- [ ] Validate species names with the Global Names Resolver database (accessed through taxize::gnr_resolve)
- [ ] Automatically correct invalid names using fuzzy matching in canonical names (challenge: avoid wrong resolution by checking that resolved taxononomic group matches expectations, they used https://www.ncbi.nlm.nih.gov/taxonomy)
- [ ] Check for synonyms in data using the ITIS database
- [ ] Manual corrections
- [ ] Afterwards, remove networks that have species in two incompatible guilds

An additional challenge that I will have is to determine the guild of each species.
For most network types, it can be easily retrieved as rows are one guild (e.g. plant) and columns another (e.g. pollinator). But I will have an issue with food webs, in which it is perfectly acceptable that a species is at the same time consumer and resource.

For now, the first step is to create a data table that contains all raw species names in the Web of Life dataset as well as some basic information on network ID, location ID and guild.

#### 16 April 2021
##### Starting to read and prepare interactions data:
I slowly started to read Web of Life raw data and think about the next steps.
Next week, I will need to start:
* Formatting data on species (as a table with species information)
* Cleaning species names (merging synonyms, resolving weird names)
* Formatting data on interactions (as a table with information on each interaction for each network)
* Format interactions metadata by adding a column for location ID and the possibility to assign locations manually

##### Jake's feedback on the outline:
Today I received feedback from Jake about the outline.
In particular, he made the following points:
* Think some more about your hypotheses regarding how the different factors I mention will influence the fidelity of interactions
* Make sure I am crystal clear about the logical flow between ideas --- the development of ideas towards the main goal of the study are not obvious yet
* Think about a plan-B in case things go wrong --- the project is quite ambitious
* Technically, in the opportunism-fidelity spectrum, the two ends are perfect non-random fidelity and perfect non-random infidelity, with random opportunism as a middle ground
* Minor point: be clear that by existence of pairwise interactions I mean the probability of an interaction between a given pair of species

##### Starting the research log:
Starting to write the research log here.
The project began on March.
Until now, I have done some literature search, read *Statistical Rethinking* --- a book on Bayesian statistics ---, written the outline, prepared the docker container, and started to download the data.
