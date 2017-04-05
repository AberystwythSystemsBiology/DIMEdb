function render_view(base_url, id) {
    var url = base_url + "api/metabolite/?id__exact=" + id;
    var x_plot = [];
    var y_plot = [];

    $.getJSON(url, function (data) {
        // Obtain metabolite dictonary
        var metabolite = data["data"][0];
        // Change the document title.
        $(document).prop('title', metabolite["name"] + " : DIMEdb - Direct Infusion MEtabolite Database");

        $("#met_name").html(metabolite["name"]);
        $("#met_name_modal").html(metabolite["name"]);
        $("#get_json").attr("href", base_url + "api/metabolite/?id__exact=" + metabolite["id"]);
        $("#molecular_formula").html(metabolite["molecular_formula"].replace(/([0-9]+)/g, '<sub>$1</sub>'));
        //$("#molecular_formula_search").attr("href", base_url + )
        $("#accurate_mass").html(metabolite["accurate_mass"].toFixed(6));
        $("#num_atoms").html(metabolite["num_atoms"]);


        $("#molecular_formula_search").click(function() {
            render_search_results("mol_form_search_results", "api/metabolites/?molecular_formula__exact=" + metabolite["molecular_formula"], 10);
            $("#molecular_formula_modal").modal("toggle");
        });

        kegg_dict = {
            "map00010": "Glycolysis / Gluconeogenesis",
            "map00020": "Citrate cycle (TCA cycle)",
            "map00030": "Pentose phosphate pathway",
            "map00040": "Pentose and glucuronate interconversions",
            "map00051": "Fructose and mannose metabolism",
            "map00052": "Galactose metabolism",
            "map00053": "Ascorbate and aldarate metabolism",
            "map00061": "Fatty acid biosynthesis",
            "map00062": "Fatty acid elongation",
            "map00071": "Fatty acid degradation",
            "map00072": "Synthesis and degradation of ketone bodies",
            "map00073": "Cutin, suberine and wax biosynthesis",
            "map00100": "Steroid biosynthesis",
            "map00120": "Primary bile acid biosynthesis",
            "map00121": "Secondary bile acid biosynthesis",
            "map00130": "Ubiquinone and other terpenoid-quinone biosynthesis",
            "map00140": "Steroid hormone biosynthesis",
            "map00190": "Oxidative phosphorylation",
            "map00195": "Photosynthesis",
            "map00196": "Photosynthesis - antenna proteins",
            "map00220": "Arginine biosynthesis",
            "map00230": "Purine metabolism",
            "map00231": "Puromycin biosynthesis",
            "map00232": "Caffeine metabolism",
            "map00240": "Pyrimidine metabolism",
            "map00250": "Alanine, aspartate and glutamate metabolism",
            "map00253": "Tetracycline biosynthesis",
            "map00254": "Aflatoxin biosynthesis",
            "map00260": "Glycine, serine and threonine metabolism",
            "map00261": "Monobactam biosynthesis",
            "map00270": "Cysteine and methionine metabolism",
            "map00280": "Valine, leucine and isoleucine degradation",
            "map00281": "Geraniol degradation",
            "map00290": "Valine, leucine and isoleucine biosynthesis",
            "map00300": "Lysine biosynthesis",
            "map00310": "Lysine degradation",
            "map00311": "Penicillin and cephalosporin biosynthesis",
            "map00330": "Arginine and proline metabolism",
            "map00331": "Clavulanic acid biosynthesis",
            "map00332": "Carbapenem biosynthesis",
            "map00340": "Histidine metabolism",
            "map00350": "Tyrosine metabolism",
            "map00351": "DDT degradation",
            "map00360": "Phenylalanine metabolism",
            "map00361": "Chlorocyclohexane and chlorobenzene degradation",
            "map00362": "Benzoate degradation",
            "map00363": "Bisphenol degradation",
            "map00364": "Fluorobenzoate degradation",
            "map00365": "Furfural degradation",
            "map00380": "Tryptophan metabolism",
            "map00400": "Phenylalanine, tyrosine and tryptophan biosynthesis",
            "map00401": "Novobiocin biosynthesis",
            "map00402": "Benzoxazinoid biosynthesis",
            "map00403": "Indole diterpene alkaloid biosynthesis",
            "map00404": "Staurosporine biosynthesis",
            "map00410": "beta-Alanine metabolism",
            "map00430": "Taurine and hypotaurine metabolism",
            "map00440": "Phosphonate and phosphinate metabolism",
            "map00450": "Selenocompound metabolism",
            "map00460": "Cyanoamino acid metabolism",
            "map00471": "D-Glutamine and D-glutamate metabolism",
            "map00472": "D-Arginine and D-ornithine metabolism",
            "map00473": "D-Alanine metabolism",
            "map00480": "Glutathione metabolism",
            "map00500": "Starch and sucrose metabolism",
            "map00510": "N-Glycan biosynthesis",
            "map00511": "Other glycan degradation",
            "map00512": "Mucin type O-glycan biosynthesis",
            "map00513": "Various types of N-glycan biosynthesis",
            "map00514": "Other types of O-glycan biosynthesis",
            "map00515": "Mannose type O-glycan biosynthesis",
            "map00520": "Amino sugar and nucleotide sugar metabolism",
            "map00521": "Streptomycin biosynthesis",
            "map00522": "Biosynthesis of 12-, 14- and 16-membered macrolides",
            "map00523": "Polyketide sugar unit biosynthesis",
            "map00524": "Neomycin, kanamycin and gentamicin biosynthesis",
            "map00525": "Acarbose and validamycin biosynthesis",
            "map00531": "Glycosaminoglycan degradation",
            "map00532": "Glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate",
            "map00533": "Glycosaminoglycan biosynthesis - keratan sulfate",
            "map00534": "Glycosaminoglycan biosynthesis - heparan sulfate / heparin",
            "map00540": "Lipopolysaccharide biosynthesis",
            "map00550": "Peptidoglycan biosynthesis",
            "map00561": "Glycerolipid metabolism",
            "map00562": "Inositol phosphate metabolism",
            "map00563": "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis",
            "map00564": "Glycerophospholipid metabolism",
            "map00565": "Ether lipid metabolism",
            "map00590": "Arachidonic acid metabolism",
            "map00591": "Linoleic acid metabolism",
            "map00592": "alpha-Linolenic acid metabolism",
            "map00600": "Sphingolipid metabolism",
            "map00601": "Glycosphingolipid biosynthesis - lacto and neolacto series",
            "map00603": "Glycosphingolipid biosynthesis - globo and isoglobo series",
            "map00604": "Glycosphingolipid biosynthesis - ganglio series",
            "map00620": "Pyruvate metabolism",
            "map00621": "Dioxin degradation",
            "map00622": "Xylene degradation",
            "map00623": "Toluene degradation",
            "map00624": "Polycyclic aromatic hydrocarbon degradation",
            "map00625": "Chloroalkane and chloroalkene degradation",
            "map00626": "Naphthalene degradation",
            "map00627": "Aminobenzoate degradation",
            "map00630": "Glyoxylate and dicarboxylate metabolism",
            "map00633": "Nitrotoluene degradation",
            "map00640": "Propanoate metabolism",
            "map00642": "Ethylbenzene degradation",
            "map00643": "Styrene degradation",
            "map00650": "Butanoate metabolism",
            "map00660": "C5-Branched dibasic acid metabolism",
            "map00670": "One carbon pool by folate",
            "map00680": "Methane metabolism",
            "map00710": "Carbon fixation in photosynthetic organisms",
            "map00720": "Carbon fixation pathways in prokaryotes",
            "map00730": "Thiamine metabolism",
            "map00740": "Riboflavin metabolism",
            "map00750": "Vitamin B6 metabolism",
            "map00760": "Nicotinate and nicotinamide metabolism",
            "map00770": "Pantothenate and CoA biosynthesis",
            "map00780": "Biotin metabolism",
            "map00785": "Lipoic acid metabolism",
            "map00790": "Folate biosynthesis",
            "map00791": "Atrazine degradation",
            "map00830": "Retinol metabolism",
            "map00860": "Porphyrin and chlorophyll metabolism",
            "map00900": "Terpenoid backbone biosynthesis",
            "map00901": "Indole alkaloid biosynthesis",
            "map00902": "Monoterpenoid biosynthesis",
            "map00903": "Limonene and pinene degradation",
            "map00904": "Diterpenoid biosynthesis",
            "map00905": "Brassinosteroid biosynthesis",
            "map00906": "Carotenoid biosynthesis",
            "map00908": "Zeatin biosynthesis",
            "map00909": "Sesquiterpenoid and triterpenoid biosynthesis",
            "map00910": "Nitrogen metabolism",
            "map00920": "Sulfur metabolism",
            "map00930": "Caprolactam degradation",
            "map00940": "Phenylpropanoid biosynthesis",
            "map00941": "Flavonoid biosynthesis",
            "map00942": "Anthocyanin biosynthesis",
            "map00943": "Isoflavonoid biosynthesis",
            "map00944": "Flavone and flavonol biosynthesis",
            "map00945": "Stilbenoid, diarylheptanoid and gingerol biosynthesis",
            "map00950": "Isoquinoline alkaloid biosynthesis",
            "map00960": "Tropane, piperidine and pyridine alkaloid biosynthesis",
            "map00965": "Betalain biosynthesis",
            "map00966": "Glucosinolate biosynthesis",
            "map00970": "Aminoacyl-tRNA biosynthesis",
            "map00980": "Metabolism of xenobiotics by cytochrome P450",
            "map00981": "Insect hormone biosynthesis",
            "map00982": "Drug metabolism - cytochrome P450",
            "map00983": "Drug metabolism - other enzymes",
            "map00984": "Steroid degradation",
            "map01010": "Overview of biosynthetic pathways",
            "map01040": "Biosynthesis of unsaturated fatty acids",
            "map01051": "Biosynthesis of ansamycins",
            "map01052": "Type I polyketide structures",
            "map01053": "Biosynthesis of siderophore group nonribosomal peptides",
            "map01054": "Nonribosomal peptide structures",
            "map01055": "Biosynthesis of vancomycin group antibiotics",
            "map01056": "Biosynthesis of type II polyketide backbone",
            "map01057": "Biosynthesis of type II polyketide products",
            "map01058": "Acridone alkaloid biosynthesis",
            "map01059": "Biosynthesis of enediyne antibiotics",
            "map01060": "Biosynthesis of plant secondary metabolites",
            "map01061": "Biosynthesis of phenylpropanoids",
            "map01062": "Biosynthesis of terpenoids and steroids",
            "map01063": "Biosynthesis of alkaloids derived from shikimate pathway",
            "map01064": "Biosynthesis of alkaloids derived from ornithine, lysine and nicotinic acid",
            "map01065": "Biosynthesis of alkaloids derived from histidine and purine",
            "map01066": "Biosynthesis of alkaloids derived from terpenoid and polyketide",
            "map01070": "Biosynthesis of plant hormones",
            "map01100": "Metabolic pathways",
            "map01110": "Biosynthesis of secondary metabolites",
            "map01120": "Microbial metabolism in diverse environments",
            "map01130": "Biosynthesis of antibiotics",
            "map01200": "Carbon metabolism",
            "map01210": "2-Oxocarboxylic acid metabolism",
            "map01212": "Fatty acid metabolism",
            "map01220": "Degradation of aromatic compounds",
            "map01230": "Biosynthesis of amino acids",
            "map01501": "beta-Lactam resistance",
            "map01502": "Vancomycin resistance",
            "map01503": "Cationic antimicrobial peptide (CAMP) resistance",
            "map01521": "EGFR tyrosine kinase inhibitor resistance",
            "map01522": "Endocrine resistance",
            "map01523": "Antifolate resistance",
            "map01524": "Platinum drug resistance",
            "map02010": "ABC transporters",
            "map02020": "Two-component system",
            "map02024": "Quorum sensing",
            "map02025": "Biofilm formation - Pseudomonas aeruginosa",
            "map02026": "Biofilm formation - Escherichia coli",
            "map02030": "Bacterial chemotaxis",
            "map02040": "Flagellar assembly",
            "map02060": "Phosphotransferase system (PTS)",
            "map03008": "Ribosome biogenesis in eukaryotes",
            "map03010": "Ribosome",
            "map03013": "RNA transport",
            "map03015": "mRNA surveillance pathway",
            "map03018": "RNA degradation",
            "map03020": "RNA polymerase",
            "map03022": "Basal transcription factors",
            "map03030": "DNA replication",
            "map03040": "Spliceosome",
            "map03050": "Proteasome",
            "map03060": "Protein export",
            "map03070": "Bacterial secretion system",
            "map03320": "PPAR signaling pathway",
            "map03410": "Base excision repair",
            "map03420": "Nucleotide excision repair",
            "map03430": "Mismatch repair",
            "map03440": "Homologous recombination",
            "map03450": "Non-homologous end-joining",
            "map03460": "Fanconi anemia pathway",
            "map04010": "MAPK signaling pathway",
            "map04011": "MAPK signaling pathway - yeast",
            "map04012": "ErbB signaling pathway",
            "map04013": "MAPK signaling pathway - fly",
            "map04014": "Ras signaling pathway",
            "map04015": "Rap1 signaling pathway",
            "map04016": "MAPK signaling pathway - plant",
            "map04020": "Calcium signaling pathway",
            "map04022": "cGMP-PKG signaling pathway",
            "map04024": "cAMP signaling pathway",
            "map04060": "Cytokine-cytokine receptor interaction",
            "map04062": "Chemokine signaling pathway",
            "map04064": "NF-kappa B signaling pathway",
            "map04066": "HIF-1 signaling pathway",
            "map04068": "FoxO signaling pathway",
            "map04070": "Phosphatidylinositol signaling system",
            "map04071": "Sphingolipid signaling pathway",
            "map04072": "Phospholipase D signaling pathway",
            "map04075": "Plant hormone signal transduction",
            "map04080": "Neuroactive ligand-receptor interaction",
            "map04110": "Cell cycle",
            "map04111": "Cell cycle - yeast",
            "map04112": "Cell cycle - Caulobacter",
            "map04113": "Meiosis - yeast",
            "map04114": "Oocyte meiosis",
            "map04115": "p53 signaling pathway",
            "map04120": "Ubiquitin mediated proteolysis",
            "map04122": "Sulfur relay system",
            "map04130": "SNARE interactions in vesicular transport",
            "map04137": "Mitophagy - animal",
            "map04138": "Autophagy - yeast",
            "map04139": "Mitophagy - yeast",
            "map04140": "Autophagy - animal",
            "map04141": "Protein processing in endoplasmic reticulum",
            "map04142": "Lysosome",
            "map04144": "Endocytosis",
            "map04145": "Phagosome",
            "map04146": "Peroxisome",
            "map04150": "mTOR signaling pathway",
            "map04151": "PI3K-Akt signaling pathway",
            "map04152": "AMPK signaling pathway",
            "map04210": "Apoptosis",
            "map04211": "Longevity regulating pathway",
            "map04212": "Longevity regulating pathway - worm",
            "map04213": "Longevity regulating pathway - multiple species",
            "map04214": "Apoptosis - fly",
            "map04215": "Apoptosis - multiple species",
            "map04260": "Cardiac muscle contraction",
            "map04261": "Adrenergic signaling in cardiomyocytes",
            "map04270": "Vascular smooth muscle contraction",
            "map04310": "Wnt signaling pathway",
            "map04320": "Dorso-ventral axis formation",
            "map04330": "Notch signaling pathway",
            "map04340": "Hedgehog signaling pathway",
            "map04341": "Hedgehog signaling pathway - fly",
            "map04350": "TGF-beta signaling pathway",
            "map04360": "Axon guidance",
            "map04370": "VEGF signaling pathway",
            "map04371": "Apelin signaling pathway",
            "map04380": "Osteoclast differentiation",
            "map04390": "Hippo signaling pathway",
            "map04391": "Hippo signaling pathway - fly",
            "map04392": "Hippo signaling pathway -multiple species",
            "map04510": "Focal adhesion",
            "map04512": "ECM-receptor interaction",
            "map04514": "Cell adhesion molecules (CAMs)",
            "map04520": "Adherens junction",
            "map04530": "Tight junction",
            "map04540": "Gap junction",
            "map04550": "Signaling pathways regulating pluripotency of stem cells",
            "map04610": "Complement and coagulation cascades",
            "map04611": "Platelet activation",
            "map04612": "Antigen processing and presentation",
            "map04614": "Renin-angiotensin system",
            "map04620": "Toll-like receptor signaling pathway",
            "map04621": "NOD-like receptor signaling pathway",
            "map04622": "RIG-I-like receptor signaling pathway",
            "map04623": "Cytosolic DNA-sensing pathway",
            "map04624": "Toll and Imd signaling pathway",
            "map04626": "Plant-pathogen interaction",
            "map04630": "Jak-STAT signaling pathway",
            "map04640": "Hematopoietic cell lineage",
            "map04650": "Natural killer cell mediated cytotoxicity",
            "map04657": "IL-17 signaling pathway",
            "map04658": "Th1 and Th2 cell differentiation",
            "map04659": "Th17 cell differentiation",
            "map04660": "T cell receptor signaling pathway",
            "map04662": "B cell receptor signaling pathway",
            "map04664": "Fc epsilon RI signaling pathway",
            "map04666": "Fc gamma R-mediated phagocytosis",
            "map04668": "TNF signaling pathway",
            "map04670": "Leukocyte transendothelial migration",
            "map04672": "Intestinal immune network for IgA production",
            "map04710": "Circadian rhythm",
            "map04711": "Circadian rhythm - fly",
            "map04712": "Circadian rhythm - plant",
            "map04713": "Circadian entrainment",
            "map04720": "Long-term potentiation",
            "map04721": "Synaptic vesicle cycle",
            "map04722": "Neurotrophin signaling pathway",
            "map04723": "Retrograde endocannabinoid signaling",
            "map04724": "Glutamatergic synapse",
            "map04725": "Cholinergic synapse",
            "map04726": "Serotonergic synapse",
            "map04727": "GABAergic synapse",
            "map04728": "Dopaminergic synapse",
            "map04730": "Long-term depression",
            "map04740": "Olfactory transduction",
            "map04742": "Taste transduction",
            "map04744": "Phototransduction",
            "map04745": "Phototransduction - fly",
            "map04750": "Inflammatory mediator regulation of TRP channels",
            "map04810": "Regulation of actin cytoskeleton",
            "map04910": "Insulin signaling pathway",
            "map04911": "Insulin secretion",
            "map04912": "GnRH signaling pathway",
            "map04913": "Ovarian steroidogenesis",
            "map04914": "Progesterone-mediated oocyte maturation",
            "map04915": "Estrogen signaling pathway",
            "map04916": "Melanogenesis",
            "map04917": "Prolactin signaling pathway",
            "map04918": "Thyroid hormone synthesis",
            "map04919": "Thyroid hormone signaling pathway",
            "map04920": "Adipocytokine signaling pathway",
            "map04921": "Oxytocin signaling pathway",
            "map04922": "Glucagon signaling pathway",
            "map04923": "Regulation of lipolysis in adipocytes",
            "map04924": "Renin secretion",
            "map04925": "Aldosterone synthesis and secretion",
            "map04930": "Type II diabetes mellitus",
            "map04931": "Insulin resistance",
            "map04932": "Non-alcoholic fatty liver disease (NAFLD)",
            "map04933": "AGE-RAGE signaling pathway in diabetic complications",
            "map04940": "Type I diabetes mellitus",
            "map04950": "Maturity onset diabetes of the young",
            "map04960": "Aldosterone-regulated sodium reabsorption",
            "map04961": "Endocrine and other factor-regulated calcium reabsorption",
            "map04962": "Vasopressin-regulated water reabsorption",
            "map04964": "Proximal tubule bicarbonate reclamation",
            "map04966": "Collecting duct acid secretion",
            "map04970": "Salivary secretion",
            "map04971": "Gastric acid secretion",
            "map04972": "Pancreatic secretion",
            "map04973": "Carbohydrate digestion and absorption",
            "map04974": "Protein digestion and absorption",
            "map04975": "Fat digestion and absorption",
            "map04976": "Bile secretion",
            "map04977": "Vitamin digestion and absorption",
            "map04978": "Mineral absorption",
            "map05010": "Alzheimer's disease",
            "map05012": "Parkinson's disease",
            "map05014": "Amyotrophic lateral sclerosis (ALS)",
            "map05016": "Huntington's disease",
            "map05020": "Prion diseases",
            "map05030": "Cocaine addiction",
            "map05031": "Amphetamine addiction",
            "map05032": "Morphine addiction",
            "map05033": "Nicotine addiction",
            "map05034": "Alcoholism",
            "map05100": "Bacterial invasion of epithelial cells",
            "map05110": "Vibrio cholerae infection",
            "map05111": "Biofilm formation - Vibrio cholerae",
            "map05120": "Epithelial cell signaling in Helicobacter pylori infection",
            "map05130": "Pathogenic Escherichia coli infection",
            "map05131": "Shigellosis",
            "map05132": "Salmonella infection",
            "map05133": "Pertussis",
            "map05134": "Legionellosis",
            "map05140": "Leishmaniasis",
            "map05142": "Chagas disease (American trypanosomiasis)",
            "map05143": "African trypanosomiasis",
            "map05144": "Malaria",
            "map05145": "Toxoplasmosis",
            "map05146": "Amoebiasis",
            "map05150": "Staphylococcus aureus infection",
            "map05152": "Tuberculosis",
            "map05160": "Hepatitis C",
            "map05161": "Hepatitis B",
            "map05162": "Measles",
            "map05164": "Influenza A",
            "map05166": "HTLV-I infection",
            "map05168": "Herpes simplex infection",
            "map05169": "Epstein-Barr virus infection",
            "map05200": "Pathways in cancer",
            "map05202": "Transcriptional misregulation in cancer",
            "map05203": "Viral carcinogenesis",
            "map05204": "Chemical carcinogenesis",
            "map05205": "Proteoglycans in cancer",
            "map05206": "MicroRNAs in cancer",
            "map05210": "Colorectal cancer",
            "map05211": "Renal cell carcinoma",
            "map05212": "Pancreatic cancer",
            "map05213": "Endometrial cancer",
            "map05214": "Glioma",
            "map05215": "Prostate cancer",
            "map05216": "Thyroid cancer",
            "map05217": "Basal cell carcinoma",
            "map05218": "Melanoma",
            "map05219": "Bladder cancer",
            "map05220": "Chronic myeloid leukemia",
            "map05221": "Acute myeloid leukemia",
            "map05222": "Small cell lung cancer",
            "map05223": "Non-small cell lung cancer",
            "map05224": "Breast cancer",
            "map05230": "Central carbon metabolism in cancer",
            "map05231": "Choline metabolism in cancer",
            "map05310": "Asthma",
            "map05320": "Autoimmune thyroid disease",
            "map05321": "Inflammatory bowel disease (IBD)",
            "map05322": "Systemic lupus erythematosus",
            "map05323": "Rheumatoid arthritis",
            "map05330": "Allograft rejection",
            "map05332": "Graft-versus-host disease",
            "map05340": "Primary immunodeficiency",
            "map05410": "Hypertrophic cardiomyopathy (HCM)",
            "map05412": "Arrhythmogenic right ventricular cardiomyopathy (ARVC)",
            "map05414": "Dilated cardiomyopathy",
            "map05416": "Viral myocarditis",
            "map05418": "Fluid shear stress and atherosclerosis",
            "map07011": "Penicillins",
            "map07012": "Cephalosporins - parenteral agents",
            "map07013": "Cephalosporins - oral agents",
            "map07014": "Quinolones",
            "map07015": "Local analgesics",
            "map07016": "Sulfonamide derivatives - sulfa drugs",
            "map07017": "Sulfonamide derivatives - diuretics",
            "map07018": "Sulfonamide derivatives - hypoglycemic agents",
            "map07019": "Tetracyclines",
            "map07020": "Macrolides and ketolides",
            "map07021": "Aminoglycosides",
            "map07023": "Rifamycins",
            "map07024": "HMG-CoA reductase inhibitors",
            "map07025": "Quinolines",
            "map07026": "Antifungal agents",
            "map07027": "Antidepressants",
            "map07028": "Antipsychotics",
            "map07029": "Phenothiazines",
            "map07030": "Anxiolytics",
            "map07031": "Butyrophenones",
            "map07032": "Hypnotics",
            "map07033": "Anticonvulsants",
            "map07034": "Eicosanoids",
            "map07035": "Prostaglandins",
            "map07036": "Calcium channel blocking drugs",
            "map07037": "Antiarrhythmic drugs",
            "map07038": "Antiulcer drugs",
            "map07039": "Opioid analgesics",
            "map07040": "Antineoplastics - alkylating agents",
            "map07041": "Antineoplastics - antimetabolic agents",
            "map07042": "Antineoplastics - agents from natural products",
            "map07043": "Antineoplastics - hormones",
            "map07044": "Antiviral agents",
            "map07045": "Antineoplastics - protein kinases inhibitors",
            "map07046": "Immunosuppressive agents",
            "map07047": "Osteoporosis drugs",
            "map07048": "Antimigraines",
            "map07049": "Antithrombosis agents",
            "map07050": "Antirheumatics - DMARDs and biological agents",
            "map07051": "Antidiabetics",
            "map07052": "Antidyslipidemic agents",
            "map07053": "Anti-HIV agents",
            "map07054": "Antiglaucoma agents",
            "map07055": "Sulfonamide derivatives - overview",
            "map07056": "Agents for Alzheimer-type dementia",
            "map07057": "Antiparkinsonian agents",
            "map07110": "Benzoic acid family",
            "map07112": "1,2-Diphenyl substitution family",
            "map07114": "Naphthalene family",
            "map07117": "Benzodiazepine family",
            "map07211": "Serotonin receptor agonists/antagonists",
            "map07212": "Histamine H1 receptor antagonists",
            "map07213": "Dopamine receptor agonists/antagonists",
            "map07214": "beta-Adrenergic receptor agonists/antagonists",
            "map07215": "alpha-Adrenergic receptor agonists/antagonists",
            "map07216": "Catecholamine transferase inhibitors",
            "map07217": "Renin-angiotensin system inhibitors",
            "map07218": "HIV protease inhibitors",
            "map07219": "Cyclooxygenase inhibitors",
            "map07220": "Cholinergic and anticholinergic drugs",
            "map07221": "Nicotinic cholinergic receptor antagonists",
            "map07222": "Peroxisome proliferator-activated receptor (PPAR) agonists",
            "map07223": "Retinoic acid receptor (RAR) and retinoid X receptor (RXR) agonists/antagonists",
            "map07224": "Opioid receptor agonists/antagonists",
            "map07225": "Glucocorticoid and meneralocorticoid receptor agonists/antagonists",
            "map07226": "Progesterone, androgen and estrogen receptor agonists/antagonists",
            "map07227": "Histamine H2/H3 receptor agonists/antagonists",
            "map07228": "Eicosanoid receptor agonists/antagonists",
            "map07229": "Angiotensin receptor and endothelin receptor antagonists",
            "map07230": "GABA-A receptor agonists/antagonists",
            "map07231": "Sodium channel blocking drugs",
            "map07232": "Potassium channel blocking and opening drugs",
            "map07233": "Ion transporter inhibitors",
            "map07234": "Neurotransmitter transporter inhibitors",
            "map07235": "N-Metyl-D-aspartic acid receptor antagonists"
        };

        if (metabolite["pathways"] == null) {

        }
        else {
            metabolite["pathways"].forEach(function (pathway, index) {
                $("#pathway_listgroup").append(
                    "<li class='list-group-item'>" + kegg_dict[pathway] +
                    "<button id='view_pathway_button' class='btn btn-success btn-sm pull-right' name='"+pathway+"_"+kegg_dict[pathway]+"'>" +
                    "<i class='glyphicon glyphicon-eye-open'></i> View</button><div class='clearfix'></div></li>"
                );
            });
        }


        for (source in metabolite["sources"]) {
            if (metabolite["sources"][source]) {
                var source_id = metabolite["sources"][source];
                if (source == "kegg_id") {
                    var kegg_url = "http://www.genome.jp/dbget-bin/www_bget?compound+" + String(source_id);
                    $("#data_sources").append("<a href='" + kegg_url + "'><button class='btn btn-success btn-sm btn-space'><i class='glyphicon glyphicon-link'></i> KEGG</button></a>");
                }

                if (source == "chebi_id") {
                    var chebi_url = "http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:" + String(source_id);
                    $("#data_sources").append("<a href='" + chebi_url + "'><button class='btn btn-danger btn-sm btn-space'><i class='glyphicon glyphicon-link'></i> CHEBI</button></a>");

                }

                if (source == "pubchem_id") {
                    var pubchem_url = "https://pubchem.ncbi.nlm.nih.gov/compound/" + String(source_id);
                    $("#data_sources").append("<a href='" + pubchem_url + "'><button class='btn btn-warning btn-sm btn-space'><i class='glyphicon glyphicon-link'></i> PubChem</button></a>");

                }

                if (source == "hmdb_id") {
                    var hmdb_url = "http://www.hmdb.ca/metabolites/" + String(source_id);
                    $("#data_sources").append("<a href='" + hmdb_url + "'><button class='btn btn-default btn-sm btn-space'><i class='glyphicon glyphicon-link'></i> HMDB</button></a>");

                }
            }
        }


        $("#smiles").html(metabolite["smiles"]);
        $("#inchi").html(metabolite["inchi"]);


        if (metabolite["synonyms"].length == 0) {
            $("#synonyms").append("<span class='text-info'>Not available</span>");

        }
        else {
            for (indx in metabolite["synonyms"]) {
                $("#synonyms").append(metabolite["synonyms"][indx] + "; ");
            }
        }


        if (metabolite["origins"] == null) {
            $("#origins").append("<i class='text-info'>Not available</i>");
        }
        else {
            for (indx in metabolite["origins"]) {
                $("#origins").append(metabolite["origins"][indx] + "; ");
            }
        }

        if (metabolite["biofluid_locations"] == null) {
            $("#biofluids").append("<span class='text-info'>Not available</span>");
        }
        else {
            for (indx in metabolite["biofluid_locations"]) {
                $("#biofluids").append(metabolite["biofluid_locations"][indx] + "; ");
            }
        }

        if (metabolite["tissue_locations"] == null) {
            $("#tissues").append("<span class='text-primary'>Not available</span>");
        }
        else {
            for (indx in metabolite["tissue_locations"]) {
                $("#tissues").append(metabolite["tissue_locations"][indx] + "; ");
            }
        }

        $.ajax({
            url: base_url + "gen_structure/" + metabolite["id"],
            success: function (result) {
                $("#structure").attr("src", "data:image/png;base64," + result);
            }
        });

        for (result in metabolite["adducts"]) {
            if (result == "neutral") {
                $("#neutral_mass").html(metabolite["adducts"]["neutral"]["peaks"][0]["accurate_mass"]);
            }

            else {
                for (adduct_idx in metabolite["adducts"][result]["peaks"]) {
                    var adduct = metabolite["adducts"][result]["peaks"][adduct_idx];
                    $("#" + result + "_adduct").append("<li class='list-group-item'><b>" + adduct["type"] + ":</b> " + adduct["accurate_mass"].toFixed(4) + "<button id='isotope' name='" + result + "_" + adduct_idx + "' class='btn btn-primary btn-sm pull-right'><i class='glyphicon glyphicon glyphicon-stats'></i> Isotope</button><div class='clearfix'></div></li>");
                }
            }
        }

        for (i in metabolite["adducts"]["neutral"]["peaks"][0]["isotopic_distribution"]) {
            var spectra = metabolite["adducts"]["neutral"]["peaks"][0]["isotopic_distribution"][i];
            x_plot.push(spectra[0]);
            y_plot.push(spectra[1]);
        }

        var isotopic_data = [{
            x: x_plot,
            y: y_plot,
            type: 'bar',
            marker: {
                color: 'rgba(0, 0, 0, 1)'
            }
        }];

        var layout = {
            xaxis: {
                title: 'Mass-to-ion (m/z)',
                showgrid: false,
                range: [Math.min.apply(Math, x_plot) - 0.5, Math.max.apply(Math, x_plot) + 1]
            },
            yaxis: {
                title: 'Relative Intensity (%)'
            },
            margin: {
                l: 50,
                r: 50,
                b: 50,
                t: 50,
                pad: 1
            },
            bargap: 0.99
        };
        Plotly.newPlot("distribution_chart", isotopic_data, layout, {displayModeBar: false});
        $("#container").show();

        $('[id="add_to_clipboard"]').click(function() {
            console.log(id);
            clipboard_hook(id);
        });

        $('[id="view_pathway_button"]').click(function() {
            var pathway_array = $(this).attr("name").split("_");
            $("#pv_name").html(pathway_array[0]+" - "+pathway_array[1]);
            var url = "http://www.genome.jp/kegg-bin/show_pathway?"+pathway_array[0]+"+"+metabolite["sources"]["kegg"];
            $("#pathway_iframe").attr("src", url);
            $("#view_pathway").modal("toggle");
        });

        $('[id="isotope"]').click(function () {
            var x_modal_plot = [];
            var y_modal_plot = [];
            var isotope_array = $(this).attr("name").split("_");
            var adduct_dict = metabolite["adducts"][isotope_array[0]]["peaks"][isotope_array[1]];

            $("#distribution_table tbody").html("");

            $("#dm_adduct").html(adduct_dict["type"]);
            for (i in adduct_dict["isotopic_distribution"]) {
                var spectra = adduct_dict["isotopic_distribution"][i];
                $("#distribution_table tbody").append("<tr><td class='text-center'>" + spectra[0].toFixed(4) + "</td><td class='text-center'>" + spectra[1].toFixed(3) + "</td></tr>");
                x_modal_plot.push(spectra[0]);
                y_modal_plot.push(spectra[1]);
            }

            var modal_isotopic_data = [{
                x: x_modal_plot,
                y: y_modal_plot,
                type: 'bar',
                marker: {
                    color: 'rgba(0, 0, 0, 1)'
                }
            }];

            var modal_layout = {
                width : 500,
                height : 350,
                xaxis: {
                    title: 'Mass-to-ion (m/z)',
                    showgrid: false,
                    range: [Math.min.apply(Math, x_modal_plot) - 0.75, Math.max.apply(Math, x_modal_plot) + 1]
                },
                yaxis: {
                    title: 'Relative Intensity (%)'
                },
                margin: {
                    l: 50,
                    r: 50,
                    b: 50,
                    t: 50,
                    pad: 1
                },
                bargap: 0.99,
                autosize: false
            };

            Plotly.newPlot("modal_distribution_chart", modal_isotopic_data, modal_layout, {displayModeBar: false});

            $("#distribution_modal").modal("toggle");
        });
    });
}