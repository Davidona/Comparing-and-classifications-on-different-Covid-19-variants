# Comparing-and-classifications-on-different-Covid-19-variants
In partial fulfillment of the requirements for Capstone Project in Software Engineering (61998)

### 4.2 Product

The outcome of the research is a software application that can find abnormalities in viruses,
which might suggest that the virus was not created in a natural manner but rather in an engineering
process. In addition to classify genomes into clusters based on similarities.

#### 4.2.1 Algorithm implementation

Two data sets of genomes are needed, one of SARS-CoV- 2 genome sequences and the
other data set is non-SARS-CoV- 2 genomes that are called imposters. To determine whether the
SARS-CoV- 2 genomes were engineeringly altered with, The SARS-CoV- 2 genome is compared
to the imposters. The comparison method between genomes will be an alignment free based
method.

##### 4.2.1.1 Distance measure algorithm

To measure the distance between two genomes (Covid and imposter) an alignment free
frequency-based method using the Canberra distance is used [see section 2.2.3].

For a string or a genome sequence G, where |G| denotes the length of G, and G(i) is the
i-th character of G, 1 ≤ i ≤ |G|. G [i...j] is the (contiguous) substring of G from i to j.

let G1...n be a SARS-CoV-2 genome in the dataset, and let G [i...j] be a genome substring of G.
there is an empty distance vector D, and G is divided into k-mers usually substrings GK, (k-mers
are substrings of length k contained within a biological sequence). GK= {G [ 0 ...k], G [k..2k], ...,
G [i...j], ..., G [|G|-k...|G|]}. Usually K-mers are substrings that act like a moving window where


k is the window size and it jumps 1 character at a time [see section 2.2.1]. In this project the

window size is still k but the window jumps 𝒌𝟐 characters at a time.

for each substring G [i...j] of size k (k-mer) in GK it is measured against all k-mers of the other
genome (imposter), using Canberra alignment free method [see section 2.2.3]. All the measures
per G [i...j] are summed and normalized by the size of the imposter and added to the distance
vector D.

##### 4.2.1.2 Clustering algorithm

for each substring in GK (where GK is list of substrings of a Covid genome) [see section
4.2.1.1] it is compared with the k-mers of both imposters. The imposters are points of reference
represented as 0 and 1 [see figure 5 ]. If substring GK(i) is closer to point of reference 0 then V(i)=
0 else V(i)= 1. The comparison method is from section 4.2.1.1.

Figure 11 : example of comparison between two imposter genomes and one SARS-CoV-2 where the substrings are of size 4 and a
sliding window of size 1. In this figure the substring of SARS-CoV-2 is CTCT.

Every G will have a list of vectors VList, after iterating through all the SARS-CoV- 2 G1...n, the
points of reference are changed with two other imposters. then the process starts over again for
each combination of imposters a new vector V is created and added to the list of vectors VList of
each genome G. Thus, each G will hold a list of vectors made out of comparisons between G and
all the possible imposters combinations (marked as |imposters combinations|).

VList1...n = {V1, V2..., V (^) |imposters combinations|}, for each SARS-CoV- 2 genome G 1 ...n.
Moving average is applied to all vectors (V) in Each VList 1 ...n, which simplifies the data by
smoothing it out and creating one flowing line [see section 2.3]. Then dynamic time warp (DTW)
algorithm [see section 2.4] is applied. To create a distance matrix DM for each imposters pair
against all Covid genomes. The DTW method used is FastDTW [ 16 ], Which is an Accurate
Dynamic Time Warping in Linear Time and Space. After that K-medoids algorithm [see section


2.5] is used on each distance matrix (DM). This will create clustering data (CD) for each pair of

imposters. CDList would be a list of all clustering data (CD [ (^0) ...|imposter combinations|]). The last step is
to combine all data collected from CDList to create one new matrix and use k-medoids on it as
well. To do that n⨯n zero matrix M is needed, where n is the number of covid genomes. This
matrix M is a counting matrix where M[i, j] is the number of times genome i was clustered with
genome j. After that to normalize M (Mnorm), each M[i, j] in M is divided with |imposters
combinations|:
𝐌𝐧𝐨𝐫𝐦=𝟏−(|𝑖𝑚𝑝𝑜𝑠𝑡𝑒𝑟𝑠 𝑐𝑜𝑚𝑏𝑖𝑛𝑎𝑡𝑖𝑜𝑛𝑠𝑀 |).
Subtracting the norm from 1 means, the lower the number the closer the genomes are to each other
(0 means always clustered together). This creates a new distance matrix of clusters Mnorm.
K-medoids is applied on Mnorm to create the last clustering results.
