# Research Log
*Alain Chauffoureaux*

#### 28 May 2021
##### First simulation attempt
The first simulation attempt went well using 43 plants, 87 pollinators and 15 sites. This is about 10 times less than the amount of plants, pollinators and sites in the dataset, so sadly it can be expected that the actual model could take up to 1000 times longer. Since this simulation took already 13290 seconds (3h46) on the cluster with 4 chains and 2000 iterations per chain, we might need to optimize it before using real data.

All 4 chains finished successfully, despite the following informational message being issued three times:
```
Chain 2 Iteration:    1 / 2000 [  0%]  (Warmup)
Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
Chain 2 Exception: Exception: binomial_logit_lpmf: Probability parameter is inf, but must be finite! (in '/tmp/RtmpNjYU18/model-62d6f6f2af5.stan', line 10, column 6 to line 14, column 8) (in '/tmp/RtmpNjYU18/model-62d6f6f2af5.stan', line 10, column 6 to line 14, column 8)
Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
```

Here is the message that was printed at the end:
```
All 4 chains finished successfully.
Mean chain execution time: 13290.0 seconds.
Total execution time: 13523.3 seconds.
109 of 4000 (3.0%) transitions hit the maximum treedepth limit of 10 or 2^10-1 leapfrog steps.
Trajectories that are prematurely terminated due to this limit will result in slow exploration.
Increasing the max_treedepth limit can avoid this at the expense of more computation.
If increasing max_treedepth does not remove warnings, try to reparameterize the model.
```

Sadly, it will be necessary to run the chain once again because I forgot to include the generated parameters in the dataset, and thus won't be able to properly investigate this chain. I need to ensure that they will be included in the next one.

#### 21 May 2021
##### Some statistics on the number of verified species:
Here are a few information I gathered about the number of species at different stages of the cleaning process:

**Web of Life raw data**
In total, there were 22838 species entries across all Web of Life networks from various functional groups.
These entries summed to a total of 14632 unique raw species names.

**Manual corrections**
In total, 566 raw species names were manually changed because their name was too long (>3 words) or they could not be automatically corrected:
- 212 species names that were undefined (could not manually be resolved)
- 354 species names that could manually be corrected (some because the name was too long)

**A priori names**
There were 14336 unique prior names after implementing manual corrections.

**Proposed names**
The prior names were further processed to propose 10238 unique taxonomic names:
- 2851 (including 212 manuals) are unidentified (e.g. raw name was 'Unidentified')
- 4930 had identified genus but unidentified species (species tag stripped)
- 34 had abbreviations that were stripped (e.g. ssp, sp, var)
These 10238 unique taxonomic names were proposed for automatic verification.

**Verified names**
Among the proposed names, 122 were found in different kingdoms, bringing the total to 10363 species.
The 10363 proposed names were verified in ITIS and GNR with the following outcome:
- 4101 were ITIS accepted
- 701 had a verified ITIS synonym
- 4262 had a perfect GNR match
- 1279 had a GNR fuzzy match
- 20 had at least a GNR fuzzy match and a manual kingdom correction (see below)

Verifying species names enabled to identify spelling mistakes and synonyms.
The 10363 proposed names were resolved as 9901 distinct taxonomic units.

**Removing taxonomic units from higher ranks**
Out of the 9901 taxonomic units, 143 had a rank which was above genus and were removed:
- The filter was applied to 7750 taxonomic units (2151 had unidentified kingdom)
- 9758 taxonomic units were left after applying this filter

The 9758 taxonomic units belonged to the following kingdoms:
|verified_kingdom|N|
|---|---|
|Plantae|2950|
|Animalia|4574|
|<NA>|2132|
|Other Eukaryota|74|
|Fungi|16|
|Bacteria|12|

**Assigning each Web of Life entry a taxonomic unit depending on functional group**
Each taxonomic name from each network (entry) was assigned a single taxonomic unit (with name and kingdom) depending on tits functional group:
- Altogether out of the 22838 entries, 19181 had a verified taxonomic unit
- Altogether, the 10238 distinct proposed names were assigned to 9660 distinct taxonomic units
- 20 species kingdoms were manually added because they were incorrect (see above)
- 43 (of the 22838) entries were removed despite matching a taxonomic unit, due to implausible kingdom (e.g. plants registered in Web of Life as animals)
- 3756 (out of 19181) verified entries had an unknown kingdom

Note that the kingdom matching/filtering was not applied to taxonomic with unknown kingdom:
- Note that 3756 (out of the 19181 verified entries) had an unknown kingdom

**Aggregating at species level**
Taxonomic units were aggregated at the species level:
- The 9660 distinct taxonomic units were aggregated into 9628 species/genus

**Removing problematic species and networks**
Problematic networks were removed:
- 11 aquatic networks (either marine or stream ecosystems) were removed
- 1 network without information about latitude/longitude was removed
- 16039 entries belonging to 9046 species/genus remained after applying this filter

At this point, no species/genus from Fungi, Bacteria, other Eukaryota and Archaea remained:
|Kingdom|N|
|---|---|
|Plantae|4251|
|Animalia|8485|
|<NA>|3303|

**Minimum number of locations**
GBIF occurrences were only downloaded for species present in more than 1 network:
- Out of the 9046 species/genus in more than 1 network, only 1905 were present in at least two locations
- Out of these 1905, 1825 could be found in GBIF

#### 20 May 2021
##### Download and cleaning GBIF occurrences:
Last week, I prepared code to download and clean GBIF occurrences data.
It should work now but requires a lot of computation time.
I need to start using the cluster to perform the more computationally intensive tasks.

##### Meeting with Bernat:
I met with Bernat today.
We discussed a lot about the project and the next steps to take.
In particular, here is a list of some things that I need to do:
* Finish preparing the data: clean GBIF occurrences and compute the climatic suitability of species in the network, compute the degree of generalism for each species
* Double check the number of species and occurrences
* Prepare sample data for prior predictive simulations

#### 12 May 2021
##### Starting to download and clean GBIF occurrences:
Now that I have cleaned all species names and prepared interactions, I can start downloading GBIF occurrence data and clean them.
The major steps will be:
* Retrieve and select GBIF keys to download
* Download GBIF occurrence data based on keys
* Extract and clean occurrence data

The goal is to finish this part before next week (or latest during next week), but that leaves very little time considering the amount of work required.

#### 29 April 2021
##### Starting manual names verification:
Today I finished the code to automatically correct taxonomic names using ITIS, GNR and NCBI.
Here are some stats about the resolution status of names:

|validity_status|N|
|---|---|
|Valid|11979|
|Unidentified|2641|
|Too long|17|

|verification_status|N|
|---|---|
|Unverified|369|
|GNR perfect match|4179|
|ITIS accepted|3929|
|GNR fuzzy match|1339|
|ITIS synonym|728|

I need to manually check the remaining 17 + 369 names that could not be resolved automatically.

#### 22 April 2021
##### Rethinking the strategy to clean species names:
I had to rethink the strategy to clean species names. Here are the main steps of the solution that I will try to implement:
* Link each Web of Life name to a proposed_name that will be used for taxonomic resolution. At the same time, detect any potential ambiguity issue. (already done)
* Create a taxonomic dictionary that will link any proposed name to a valid taxon (one that is present and valid in the ITIS database). This step will imply to resolve names that are ITIS synonyms, autocorrection of misspelled names using GNR, as well as some manual corrections whenever applicable (or alternatively they can be done forthe proposed name directly, the solution I have now).
* Create a database that gives information on all valid taxa from the taxonomic dictionary. This includes the rank of the taxon, the name of its higher rank, species and subspecies (when given), as well as links to parent taxonomic units.
* Use the taxonomic dictionary to give a taxonomic ID to all Web of Life names. Additionally, give a wol_species_ID and a wol_subspecies_id that will enable to distinguish ambiguous (different species/subspecies that belong to the same higher-level taxa) and unverified rows (entities that do not have matching taxonomic information)

I think that there could be a few different advantages to this new strategy:
* It makes a clear distinction between taxonomic information (that is well resolved) and entities (that can be unresolved or be different although linked to the exact same taxonomic information).
* Later on, I can use all synonyms (i.e. names that resolve to the same taxon) of a given taxon when querying other databases such as GBIF. That way, data should be found even if the database uses unusual names.
* I can decide later on at which level I want to collapse information. I can keep subspecies resolution, aggregate subspecies at species level, or even aggregate species at genus (or higher) level.

Let's hope this strategy will prove effective, but it is probable that it changes again while I implement it.

#### 21 April 2021
##### What to do with undefined species of same genus in a network:
I am not sure yet if and how I will include species that are only resolved at the genus level.
There are many, so it could make sense to still include observations of defined genus but undefined species (handling genus as if it was a species).
For now, I will also resolve genera names to leave the possibility open, and I can decide later.
That mean that it is possible that:
* In the same network, several unknown species of the same genus are both resolved to the same genus (I added the *needs_distinction* flag to the species dictionnary when it is the case)
* In the end, some species in the a genus could be defined at species level while others resolve to the genus name (I can check later if it is really the case)

#### 20 April 2021
##### To do next:
Today, I have started to prepare the species names for checking. Before performing the validation using GNR and ITIS, I need to finish the preparation step by:
* Flag unidentified species as such
* Resolve genus, species, subspecies in proposed names

##### Strategy to clean species names:
After looking at F. Cagua's code, I have thought of the strategy that I can adopt to clean species names. It is likely to change depending on the new challenges I face while I implement it. For now, I am considering to do it as follows:
* Starting point: data.table with raw species names (as used in Web of Life), network name, interaction type, location ID, functional group
* Apply manual name corrections
* Remove abbreviations
* Analyse raw name to determine genus, species and subspecies and indicate rank
* Flag names that indicate unidentified species
* Query cache to resolve names in cache without performing the next steps
* Check names that are not in the cache, proceed as follows:
  * Find out whether the name is in the ITIS database
  * Try to query GNR to resolve any spelling mistake if the name was not found in ITIS (verify that the kingdom is expected)
  * Determine the corresponding valid taxon based on ITIS synonyms
* Output: data.table with all raw species names either resolved or unresolved (with reason indicated)

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
