# TBC2026 POSTER SCRIPT
## "A, T, G, GPT: Teaching Machines the Language of Life"
### Aayan Alwani | Grade __ | __________ High School
### Category 3: Platform Tools & Technologies — Computational Biology

---

## HOOK (under the title)

In 2003, it took 13 years and $2.7 billion to read one human genome. Today, AI models trained on 9.3 trillion DNA letters from every domain of life are learning to read genomes the way ChatGPT reads English — and they're beginning to write new ones.

---

## I. TOPIC BACKGROUND

### What Are Genomic Foundation Models?

Genomic foundation models are large artificial intelligence systems trained on billions of DNA sequences from thousands of species across all domains of life — bacteria, plants, animals, and humans. Just as ChatGPT learned the rules of English by reading massive amounts of text from the internet, these models learn the "grammar" of genomes: which patterns regulate genes, how mutations affect function, and what makes the DNA of a human different from that of a bacterium.

### How Do They Work?

The core idea is called next-token prediction. Given a stretch of DNA — a sequence of the four nucleotide letters A, T, G, and C — the model predicts what nucleotide comes next. By training on trillions of bases from genomes spanning all of life on Earth, the AI builds a deep understanding of how DNA works, without ever being explicitly taught biology. It discovers the rules on its own, just by reading enough examples.

This is the same principle behind ChatGPT: predict the next word in a sentence, and in doing so, learn language. Genomic foundation models predict the next nucleotide in a genome, and in doing so, learn the language of life.

### Why Does This Matter?

These models can predict whether a genetic mutation causes disease, design synthetic genes with specific functions, annotate newly sequenced genomes, and transfer biological knowledge across species — all without running a single wet-lab experiment. While every AI prediction still requires experimental validation, the speed of biological discovery is being transformed. Tasks that once took years of laboratory work can now be narrowed down in hours.

### ChatGPT vs. Evo 2 — Side by Side

| ChatGPT | Evo 2 |
|---------|-------|
| Trained on text from the internet | Trained on DNA from all life on Earth |
| Learns grammar and meaning | Learns genetic "grammar" and regulation |
| Predicts the next word | Predicts the next nucleotide |
| Generates new text | Generates new DNA sequences |

---

## II. TECHNOLOGY TIMELINE

**1990 — Human Genome Project Launches**
A 13-year, $2.7 billion international effort to sequence the entire human genome begins, led by James Watson and Francis Collins at the National Institutes of Health.

**2003 — Human Genome Completed**
The Human Genome Project is declared complete, producing the first full sequence of human DNA. This milestone ushers in the era of genomics, but scientists quickly realize that reading DNA is far easier than understanding it.

**2005 — Next-Generation Sequencing Revolution**
Illumina and other companies introduce next-generation sequencing technology, slashing the cost and time of genome sequencing by over 1,000x. Sequencing a genome drops from billions of dollars to thousands, generating an explosion of genomic data.

**2017 — The Transformer Architecture Is Invented**
Researchers at Google publish "Attention Is All You Need," introducing the Transformer — the neural network architecture that would power ChatGPT, DALL-E, and eventually, genomic AI. This single invention becomes the foundation of the modern AI revolution.

**2021 — DNABERT: The First DNA Language Model**
Zhihan Zhou and colleagues at Mila (Montreal) publish DNABERT, the first BERT-style language model applied to DNA sequences. It proves that natural language processing techniques transfer directly to genomics — DNA can be treated as a language.

**2023 — Nucleotide Transformer Scales Up**
InstaDeep (acquired by BioNTech for ~$680 million) releases the Nucleotide Transformer, scaling DNA language models to 2.5 billion parameters and training on genomes from over 3,200 diverse humans and 850 species. It is deployed commercially in mRNA vaccine design.

**2024 — Evo: A Breakthrough in Genomic AI**
Patrick Hsu, Brian Hie, and Eric Nguyen at the Arc Institute publish Evo in the journal Science. With 7 billion parameters trained on 2.7 million prokaryotic genomes, Evo can predict gene essentiality, generate functional CRISPR systems, and design realistic phage genomes — all from a single model.

**2025 — Evo 2: The Largest Genomic Model Ever Built**
The Arc Institute releases Evo 2, scaling to 40 billion parameters trained on 9.3 trillion nucleotides from all domains of life. With a context window of 1 million base pairs, it can process entire gene clusters and predict variant pathogenicity across the human genome with zero-shot accuracy. It is open-sourced for researchers worldwide.

---

## III. BIOTECH INNOVATORS & ECONOMIC IMPACT

### Innovator 1: Evo & Evo 2 — Arc Institute

The Arc Institute, a non-profit research organization co-founded by Stripe CEO Patrick Collison, Stanford biologist Patrick Hsu, and biochemist Silvana Konermann, was established in 2022 with a $650 million endowment — one of the largest private research endowments in recent history. Core investigators Patrick Hsu and Brian Hie led the development of Evo (2024) and Evo 2 (2025), the most powerful genomic foundation models ever built. Critically, both models are fully open-source, allowing researchers worldwide to use them at no cost. Evo 2 can predict the effects of genetic mutations, design novel gene sequences, and understand genome organization across bacteria, plants, and humans. By making this technology freely available, the Arc Institute is democratizing access to cutting-edge genomic AI.

### Innovator 2: Nucleotide Transformer — BioNTech / InstaDeep

InstaDeep, an AI company founded in Tunis and London, developed the Nucleotide Transformer — a series of DNA language models scaling up to 2.5 billion parameters. The models excel at predicting gene regulatory elements, splice sites, and variant effects. In 2023, German pharmaceutical giant BioNTech acquired InstaDeep for approximately $680 million, one of the largest acquisitions in the AI-biotech space. The Nucleotide Transformer is now integrated into BioNTech's pipeline for designing mRNA vaccines and personalized cancer therapeutics, making it one of the most commercially deployed genomic foundation models in the world.

### Innovator 3: BioNeMo — NVIDIA

NVIDIA's BioNeMo platform provides the GPU infrastructure that powers genomic AI training and deployment across the pharmaceutical industry. Rather than building its own genomic model, NVIDIA enables the entire ecosystem — hosting models like the Nucleotide Transformer on optimized hardware and providing inference microservices for drug discovery companies. Partners include Genentech/Roche, Amgen, and Recursion Pharmaceuticals. The broader AI-in-genomics market is estimated at $1.5–2 billion in 2024 and projected to reach $5–10 billion by 2030, driven by demand for faster drug discovery, improved genetic diagnostics, and synthetic biology applications.

---

## IV. ETHICAL, LEGAL & SOCIAL ISSUES

### Biosecurity: The Dual-Use Dilemma

Genomic foundation models that can generate realistic DNA sequences raise serious biosecurity concerns. A model capable of designing functional genes could, in theory, be misused to engineer novel pathogens or enhance existing ones. The 2024 U.S. Executive Order on AI specifically identified risks from generative biology as a national security priority. The Arc Institute addressed these concerns by conducting red-teaming exercises, consulting with biosecurity experts, and implementing responsible-release policies before publishing Evo and Evo 2. DNA synthesis companies also maintain screening protocols through the International Gene Synthesis Consortium, but the gap between what AI can design and what screening can catch continues to narrow.

### Genomic Bias: Whose DNA Counts?

Training data for genomic foundation models is heavily skewed toward genomes of European ancestry. This means models may produce less accurate predictions for people of African, Asian, Indigenous, and other underrepresented ancestries — the very populations that already face the greatest health disparities. If these models are used in clinical genomics to diagnose disease or guide treatment, biased predictions could amplify existing inequities in healthcare. Addressing this requires intentional diversification of training datasets and equitable representation in biobank collections.

### Privacy and Consent: Outdated Frameworks

The genomic data used to train these models often comes from biobanks — large repositories of biological samples collected from patients and research volunteers. Many of these samples were collected under consent frameworks that predate the era of AI. Participants agreed to have their DNA used for research, but they could not have anticipated that their genetic information would be fed into a machine learning model capable of generating new sequences. Additionally, genomic data is inherently re-identifiable — even "anonymized" genetic data can be linked back to individuals. The Genetic Information Nondiscrimination Act (GINA) provides some protection against discrimination by employers and health insurers, but significant gaps remain: GINA does not cover life insurance, long-term care insurance, or disability insurance.

### Intellectual Property: Who Owns AI-Generated DNA?

If an AI model designs a novel gene sequence with therapeutic potential, who holds the patent — the developers of the model, the researchers who prompted it, or no one? Current patent law in the United States generally requires a human inventor, creating uncertainty around AI-generated biological materials. The legal landscape is further complicated by the fact that many foundational models (like Evo 2) are open-source, meaning the "inventions" they produce are built on publicly shared tools. Regulatory frameworks in both the U.S. and the European Union remain fragmented and largely untested in the context of generative biology.

---

## V. REFERENCES

1. Nguyen, E., Poli, M., Durrant, M.G., Thomas, B., Kang, K., Sullivan, J., ... & Hie, B., Hsu, P.D. (2024). "Sequence modeling and design from molecular to genome scale with Evo." *Science*, 386(6723).

2. Nguyen, E., Poli, M., Hie, B., Hsu, P.D. et al. (2025). "Genome modeling and design across all domains of life with Evo 2." *bioRxiv* preprint.

3. Dalla-Torre, H., Gonzalez, L., Mendoza-Revilla, J. et al. (2023). "The Nucleotide Transformer: Building and evaluating robust foundation models for human genomics." *Nature Machine Intelligence*, 5.

4. Zhou, Z., Ji, Y., Li, W. et al. (2023). "DNABERT-2: Efficient foundation model and benchmark for multi-species genomes." *arXiv* preprint.

5. National Human Genome Research Institute. (2025). "The Human Genome Project." genome.gov.

Generative AI tools (Claude, Anthropic) were used to assist with research and content development. All facts were independently verified against primary sources.
