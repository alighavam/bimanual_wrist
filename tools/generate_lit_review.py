#!/usr/bin/env python3
"""
Generate a comprehensive literature review PDF on bimanual coordination
and ipsilateral/contralateral motor representations.
Papers ranked by relevance to the bimanual wrist manuscript.
"""

from fpdf import FPDF


class LitReviewPDF(FPDF):
    @staticmethod
    def _ascii_safe(text: str) -> str:
        # fpdf core fonts are not Unicode; normalize common characters we use.
        return (
            text.replace("\u2013", "-")  # en dash
                .replace("\u2014", "-")  # em dash
                .replace("\u2212", "-")  # minus sign
                .replace("\u00d7", "x")  # multiplication sign
                .replace("\u00b0", "deg")  # degree sign
                .replace("\u00a0", " ")  # non-breaking space
        )

    def header(self):
        if self.page_no() == 1:
            return
        self.set_font("Helvetica", "I", 8)
        self.set_text_color(120, 120, 120)
        self.cell(
            0, 8,
            "Literature Review: Bimanual Coordination & Motor Representations",
            align="C",
        )
        self.ln(10)

    def footer(self):
        self.set_y(-15)
        self.set_font("Helvetica", "I", 8)
        self.set_text_color(120, 120, 120)
        self.cell(0, 10, f"Page {self.page_no()}/{{nb}}", align="C")

    def add_title_page(self):
        self.add_page()
        self.ln(40)
        self.set_font("Helvetica", "B", 22)
        self.set_text_color(30, 30, 100)
        self.cell(0, 12, "Literature Review", align="C",
                  new_x="LMARGIN", new_y="NEXT")
        self.ln(5)
        self.set_font("Helvetica", "", 14)
        self.set_text_color(60, 60, 60)
        self.cell(0, 8, "Bimanual Coordination, Ipsilateral Motor",
                  align="C", new_x="LMARGIN", new_y="NEXT")
        self.cell(0, 8, "Representations, and Reference Frames",
                  align="C", new_x="LMARGIN", new_y="NEXT")
        self.cell(0, 8, "in Motor Cortex",
                  align="C", new_x="LMARGIN", new_y="NEXT")
        self.ln(10)
        self.set_font("Helvetica", "I", 11)
        self.set_text_color(100, 100, 100)
        self.cell(0, 7, "Compiled for: Bimanual Wrist Reaching Manuscript",
                  align="C", new_x="LMARGIN", new_y="NEXT")
        self.cell(0, 7, "April 2026",
                  align="C", new_x="LMARGIN", new_y="NEXT")
        self.ln(10)
        self.set_font("Helvetica", "", 10)
        self.set_text_color(80, 80, 80)
        self.multi_cell(
            0, 6,
            "This review covers 45 papers spanning monkey electrophysiology, "
            "human fMRI, behavioral studies, and patient work on bimanual "
            "coordination, intrinsic vs. extrinsic reference frames, "
            "ipsilateral motor representations, and multivariate encoding "
            "of hand movements. Papers are ranked by relevance to the "
            "manuscript.",
            align="C",
        )

    def add_section_header(self, text):
        self.ln(4)
        self.set_font("Helvetica", "B", 12)
        self.set_text_color(30, 30, 100)
        self.multi_cell(0, 7, text)
        self.set_draw_color(30, 30, 100)
        self.line(
            self.l_margin, self.get_y(),
            self.w - self.r_margin, self.get_y(),
        )
        self.ln(4)

    def add_paper(self, number, title, authors, year, journal, doi_or_link,
                  gap, task, results, relation_score=None):
        # Check space
        if self.get_y() > 215:
            self.add_page()

        # Number + title
        self.set_font("Helvetica", "B", 10)
        self.set_text_color(30, 30, 100)
        self.multi_cell(0, 5.5, self._ascii_safe(f"{number}. {title}"))
        self.ln(1)

        # Authors / year / journal (+ relation score)
        self.set_font("Helvetica", "", 8)
        self.set_text_color(80, 80, 80)
        meta = f"{authors} ({year}). {journal}."
        if relation_score is not None:
            meta += f" Relevance: {relation_score}/10."
        self.multi_cell(0, 4.5, self._ascii_safe(meta))
        self.ln(1)

        # Link
        self.set_font("Helvetica", "U", 7.5)
        self.set_text_color(0, 0, 200)
        self.multi_cell(0, 4.5, self._ascii_safe(doi_or_link))
        self.ln(2)

        # Gap
        self.set_font("Helvetica", "", 8.5)
        self.set_text_color(40, 40, 40)
        self.multi_cell(0, 4.5, self._ascii_safe(f"[Gap] {gap}"))
        self.ln(1)

        # Task
        self.set_font("Helvetica", "", 8.5)
        self.set_text_color(40, 40, 40)
        self.multi_cell(0, 4.5, self._ascii_safe(f"[Task & Modality] {task}"))
        self.ln(1)

        # Results
        self.set_font("Helvetica", "", 8.5)
        self.set_text_color(40, 40, 40)
        self.multi_cell(0, 4.5, self._ascii_safe(f"[Key Results] {results}"))
        self.ln(2)

        # Separator
        self.set_draw_color(200, 200, 200)
        y = self.get_y()
        self.line(self.l_margin + 30, y, self.w - self.r_margin - 30, y)
        self.ln(4)


# ------------------------------------------------------------------ #
#  Paper data
# ------------------------------------------------------------------ #

def get_papers():
    return [
        # ===== SECTION 1: Bimanual Activity Decomposition & Neural Encoding =====
        {
            "title": "Neural organization of hierarchical motor sequence representations in the human neocortex",
            "authors": "Yokoi A, Diedrichsen J",
            "year": "2019",
            "journal": "Neuron, 103(6), 1178-1190.e7",
            "doi_or_link": "https://doi.org/10.1016/j.neuron.2019.06.017",
            "relation_score": 2,
            "gap": (
                "Although it is widely accepted that movement sequences are represented "
                "hierarchically (individual movements, chunks, and whole sequences), it "
                "was unclear how these different levels map onto human neocortex. In "
                "particular, do different cortical regions represent distinct levels of "
                "the hierarchy (anatomical separation), or do chunk and sequence "
                "representations coexist within the same regions?"
            ),
            "task": (
                "Participants learned and executed 8 trained right-hand finger-press "
                "sequences explicitly organized into four chunks (2–3 presses each). "
                "Used representational fMRI analyses (multivariate pattern similarity) "
                "to test for representations of individual finger presses, chunks, and "
                "whole sequences across cortex."
            ),
            "results": (
                "Found clear evidence for all three representational levels. "
                "Anatomically, individual movements were uniquely represented in "
                "primary motor cortex (M1), whereas movement chunks and entire "
                "sequences were jointly represented with substantial overlap in "
                "premotor and parietal cortices. These findings challenge a strict "
                "anatomical hierarchy where each level is cleanly separated into "
                "different regions, and instead highlight a special distinction "
                "between representations of individual movements and sequential context."
            ),
        },
        {
            "title": "Effector-invariant movement encoding in the human motor system",
            "authors": "Haar S, Dinstein I, Shelef I, Donchin O",
            "year": "2017",
            "journal": "Journal of Neuroscience, 37(37), 9054-9063",
            "doi_or_link": "https://doi.org/10.1523/JNEUROSCI.1663-17.2017",
            "relation_score": 9,
            "gap": (
                "Whether ipsilateral motor regions encode movements in the same "
                "coordinate system as contralateral movements (extrinsic/spatial "
                "vs. intrinsic/joint-based) was unclear. Understanding this "
                "relationship reveals how the brain manages bilateral motor control "
                "without causing interference."
            ),
            "task": (
                "Center-out reaching with each arm to 4 targets, measured with fMRI "
                "and multivariate pattern analysis (crossnobis dissimilarity)."
            ),
            "results": (
                "Ipsilateral and contralateral movements with symmetric joint "
                "configurations were encoded similarly by neural populations in M1, "
                "PMd, SMA, and SPL, providing evidence for effector-invariant "
                "encoding in intrinsic/joint coordinates. This is directly relevant "
                "to the manuscript's finding that contra- and ipsilateral activity "
                "patterns are intrinsically correlated rather than extrinsically "
                "correlated. The study used the same crossnobis metric as the manuscript."
            ),
        },
        {
            "title": "Future movement plans interact in sequential arm movements",
            "authors": "Kashefi M, Reschechtko S, Ariani G, Shahbazi M, Tan A, Diedrichsen J, Pruszynski JA",
            "year": "2024",
            "journal": "eLife, 13, e94485",
            "doi_or_link": "https://doi.org/10.7554/eLife.94485",
            "relation_score": 2,
            "gap": (
                "Real-world actions require planning future movements while executing "
                "the current one. It was unclear (1) how many future reaches can be "
                "planned simultaneously in longer sequences and (2) whether planning "
                "processes for multiple future movements are independent or interact "
                "with each other and with ongoing control."
            ),
            "task": (
                "Behavioral continuous sequential reaching: participants executed "
                "14-reach sequences in a planar robotic exoskeleton while the number "
                "of visible future targets (planning horizon H1–H5) and dwell times "
                "(75/200/400 ms) were manipulated. Used target-jump perturbations to "
                "probe whether +2 reaches are planned and whether future plans interact."
            ),
            "results": (
                "Participants planned at least two future reaches during execution of "
                "the current reach. Jumping the +2 target revealed partial commitment "
                "to the pre-jump +2 location, indicating advance planning of +2. "
                "Critically, planning processes interacted: correcting an ongoing reach "
                "to a +1 target jump was slower/longer when participants could also plan "
                "+2 (H3) compared to when only +1 could be planned (H2). Additionally, "
                "current reach curvature was influenced by the next reach only when the "
                "planning processes temporally overlapped, consistent with resource "
                "sharing or interaction between future movement plans."
            ),
        },
        {
            "title": "Where one hand meets the other: limb-specific and action-dependent movement plans decoded from preparatory signals in single human frontoparietal brain areas",
            "authors": "Gallivan JP, McLean DA, Flanagan JR, Culham JC",
            "year": "2013",
            "journal": "Journal of Neuroscience, 33(5), 1991-2008",
            "doi_or_link": "https://doi.org/10.1523/JNEUROSCI.0541-12.2013",
            "relation_score": 5,
            "gap": (
                "While contralateral limb planning is well-established in "
                "frontoparietal cortex, the degree to which ipsilateral limb actions "
                "are represented was underappreciated. Can fMRI multivariate methods "
                "decode which hand will be used and what action will be performed "
                "from preparatory activity?"
            ),
            "task": (
                "Reaching vs. grasping with left or right hand while maintaining "
                "fixation. Decoded from fMRI preparatory signals using MVPA."
            ),
            "results": (
                "Found much stronger ipsilateral limb representation than expected. "
                "Subregions of posterior parietal cortex, PMd, PMv, pre-SMA, and "
                "even M1 could decode upcoming ipsilateral hand actions. This "
                "challenges the strict contralateral view and supports the "
                "manuscript's finding that ipsilateral encoding is present even "
                "in M1 and S1. Premotor and parietal areas showed approximately "
                "symmetric bilateral representations."
            ),
        },
        {
            "title": "Skill learning strengthens cortical representations of motor sequences",
            "authors": "Wiestler T, Diedrichsen J",
            "year": "2013",
            "journal": "eLife, 2, e00801",
            "doi_or_link": "https://doi.org/10.7554/eLife.00801",
            "relation_score": 3,
            "gap": (
                "Whether multivariate fMRI methods can reveal fine-grained movement "
                "representations in motor cortex, beyond what univariate activation "
                "shows, and whether learning strengthens these representations needed "
                "validation."
            ),
            "task": (
                "Finger sequence learning with fMRI. Participants trained for 4 days "
                "on four left-hand finger sequences and then performed trained vs. "
                "untrained sequences during scanning. Used multivariate pattern "
                "analysis (MVPA) to discriminate sequence-specific activity patterns "
                "and quantify how reliably sequences could be distinguished."
            ),
            "results": (
                "Motor-skill acquisition was associated with the emergence of more "
                "distinguishable sequence-specific activity patterns for trained "
                "sequences, even without increases in spatially averaged activity. "
                "Both trained and untrained sequences could be discriminated in "
                "primary and secondary motor areas, but trained sequences were "
                "classified more reliably (especially in SMA). This demonstrates "
                "that learning can strengthen representational specificity even when "
                "univariate activation does not increase."
            ),
        },
        {
            "title": "Distinct representation of ipsilateral hand movements in sensorimotor areas",
            "authors": "Bruurmijn MLCM, Raemaekers M, Branco MP, Ramsey NF, Vansteensel MJ",
            "year": "2021",
            "journal": "European Journal of Neuroscience, 54(10), 7599-7608",
            "doi_or_link": "https://doi.org/10.1111/ejn.15501",
            "relation_score": 6,
            "gap": (
                "While contralateral motor cortex involvement is well-established, "
                "the role of ipsilateral sensorimotor cortex during unilateral "
                "movements remained debated. Can different ipsilateral hand gestures "
                "be decoded, and from which sub-regions of sensorimotor cortex?"
            ),
            "task": (
                "Multiple hand gestures during 7T fMRI. Multivariate "
                "decoding of ipsilateral vs. contralateral hand movement identity."
            ),
            "results": (
                "Ipsilateral hand gestures could be distinguished from each other "
                "using multivariate pattern analysis. Ipsilateral activation was "
                "strongest in anterior precentral gyrus and posterior postcentral "
                "gyrus. Critically, ipsilateral patterns were distinct from "
                "contralateral patterns, suggesting the same cortical region uses "
                "different population codes for each hand. Relevant to the "
                "manuscript's finding of ipsilateral encoding across motor regions."
            ),
        },
        {
            "title": "The role of human primary motor cortex in the production of skilled finger sequences",
            "authors": "Yokoi A, Arbuckle SA, Diedrichsen J",
            "year": "2018",
            "journal": "Journal of Neuroscience, 38(6), 1430-1442",
            "doi_or_link": "https://doi.org/10.1523/JNEUROSCI.2798-17.2017",
            "relation_score": 4,
            "gap": (
                "Although M1 is essential for producing individuated finger movements, "
                "it was unclear whether M1 also represents learned sequences of multiple "
                "finger movements (i.e., sequential context), or whether sequence "
                "representations reside primarily in premotor/parietal areas with M1 "
                "reflecting only execution of component finger presses."
            ),
            "task": (
                "Participants practiced finger press sequences and performed them during "
                "3T and 7T fMRI. Used multivariate representational analyses to test for "
                "sequence representations in M1 versus premotor/parietal cortex. "
                "Compared multi-finger sequence patterns to patterns for constituent "
                "single-finger movements, and used passive replay to test whether "
                "sequence-pattern effects reflect execution rather than hemodynamic "
                "nonlinearities."
            ),
            "results": (
                "After intensive practice, premotor and parietal areas encoded the "
                "different movement sequences, but there was little or no evidence for "
                "sequence representations in M1. Instead, M1 activity patterns during "
                "sequences were fully explained by a linear combination of patterns for "
                "individual finger movements, with a disproportionately strong weight "
                "on the first finger of the sequence (the 'first-finger effect'). Passive "
                "replay showed this effect is linked to active execution processes rather "
                "than fMRI hemodynamic nonlinearity. The results suggest M1 primarily "
                "reflects execution of component presses, while sequence context is "
                "represented in premotor/parietal cortex."
            ),
        },
        {
            "title": "Two distinct ipsilateral cortical representations for individuated finger movements",
            "authors": "Diedrichsen J, Wiestler T, Krakauer JW",
            "year": "2013",
            "journal": "Cerebral Cortex, 23(6), 1362-1377",
            "doi_or_link": "https://doi.org/10.1093/cercor/bhs120",
            "relation_score": 10,
            "gap": (
                "Ipsilateral motor cortex shows activity changes during unimanual "
                "movements, but its functional significance was unclear. Specifically, "
                "do ipsilateral cortical areas carry finger-specific information, and "
                "does the nature of ipsilateral representations change between "
                "unimanual and bimanual contexts (i.e., when both hemispheres are "
                "actively engaged)?"
            ),
            "task": (
                "Human fMRI with multivoxel pattern analysis of fine-grained activation "
                "patterns during paced isometric finger presses. Experiment 1 tested "
                "unimanual left vs. right finger presses (all 10 digits). Experiment 2 "
                "tested unimanual and bimanual presses (digits 1, 3, 5) to examine how "
                "contra- and ipsilateral finger representations interact during bimanual "
                "actions."
            ),
            "results": (
                "Showed two fundamentally different ipsilateral representations. "
                "During unimanual ipsilateral presses, primary sensory and motor "
                "cortices contained finger-specific patterns that were nearly identical "
                "to contralateral mirror-symmetric patterns, despite overall suppression; "
                "this mirrored component vanished during bimanual actions. A second "
                "ipsilateral representation emerged during bimanual actions in caudal "
                "premotor and anterior parietal cortices, where ipsilateral actions were "
                "encoded as a nonlinear modulation of patterns related to contralateral "
                "actions. Together, the results show that ipsilateral representations "
                "change their informational content and likely functional role depending "
                "on behavioral context."
            ),
        },
        {
            "title": "Hand use predicts the structure of representations in sensorimotor cortex",
            "authors": "Ejaz N, Hamada M, Diedrichsen J",
            "year": "2015",
            "journal": "Nature Neuroscience, 18(7), 1034-1040",
            "doi_or_link": "https://doi.org/10.1038/nn.4038",
            "relation_score": 5,
            "gap": (
                "What determines the structure of movement representations in "
                "sensorimotor cortex? Is it anatomy, biomechanics, or the "
                "statistics of natural hand use? Understanding this is essential "
                "for interpreting the intrinsic organization found between hands."
            ),
            "task": (
                "Single-finger presses measured with fMRI. RSA comparing brain "
                "patterns to models based on muscle anatomy and natural hand-use "
                "statistics."
            ),
            "results": (
                "The structure of finger representations in M1 was best predicted "
                "by the statistics of natural hand use (how often fingers are used "
                "together in daily life), not by muscle anatomy alone. This implies "
                "that intrinsic organization may reflect learned co-activation "
                "patterns. Relevant to interpreting why the manuscript finds "
                "intrinsic rather than extrinsic correlations between hands."
            ),
        },
        {
            "title": "Representational models: a common framework for understanding encoding, pattern-component, and representational-similarity analysis",
            "authors": "Diedrichsen J, Kriegeskorte N",
            "year": "2017",
            "journal": "PLoS Computational Biology, 13(4), e1005508",
            "doi_or_link": "https://doi.org/10.1371/journal.pcbi.1005508",
            "relation_score": 6,
            "gap": (
                "Multiple multivariate fMRI analysis methods (encoding models, RSA, "
                "pattern component models) existed but their relationships were "
                "unclear. A unifying framework was needed to understand when these "
                "methods agree or diverge."
            ),
            "task": (
                "Theoretical/methodological paper with simulations and example "
                "fMRI data. Provides the analytical framework for decomposing "
                "activity patterns into components."
            ),
            "results": (
                "Encoding models, RSA, and pattern component models (PCM) are "
                "mathematically related and can be understood as different views of "
                "the same underlying second-moment matrix of activity patterns. "
                "This framework underpins the crossnobis dissimilarity and "
                "decomposition approach used in the manuscript. Shows how to "
                "properly test whether bimanual patterns can be decomposed into "
                "additive components."
            ),
        },

        # ===== SECTION 2: Intrinsic vs. Extrinsic Reference Frames =====
        {
            "title": "Exploring interlimb constraints during bimanual graphic performance: effects of muscle grouping and direction",
            "authors": "Swinnen SP, Jardin K, Verschueren S, Meulenbroek R, Franz L, Dounskaia N, Walter CB",
            "year": "1998",
            "journal": "Behavioural Brain Research, 90(1), 79-87",
            "doi_or_link": "https://doi.org/10.1016/S0166-4328(97)00083-1",
            "relation_score": 9,
            "gap": (
                "When people perform bimanual movements, interference can arise "
                "from similarity in muscle activation patterns (intrinsic) or "
                "similarity in movement direction (extrinsic). The relative "
                "contributions of these two constraint systems to bimanual "
                "interference were debated."
            ),
            "task": (
                "Bimanual graphic (line-drawing) tasks where intrinsic "
                "(muscle) and extrinsic (spatial) congruency were independently "
                "manipulated. Behavioral kinematics."
            ),
            "results": (
                "Both intrinsic and extrinsic constraints affect bimanual "
                "coordination, but their relative dominance depends on task context. "
                "For rhythmic bimanual movements, intrinsic constraints (homologous "
                "muscle activation) were generally dominant, while extrinsic "
                "constraints became more important with visual feedback. This work "
                "established the intrinsic-extrinsic framework that the manuscript "
                "uses to test how contra- and ipsilateral activity patterns relate."
            ),
        },
        {
            "title": "On the coordination of two-handed movements",
            "authors": "Kelso JAS, Southard DL, Goodman D",
            "year": "1979",
            "journal": "Journal of Experimental Psychology: Human Perception and Performance, 5(2), 229-238",
            "doi_or_link": "https://doi.org/10.1037/0096-1523.5.2.229",
            "relation_score": 6,
            "gap": (
                "The fundamental principles governing bimanual coordination "
                "stability were not understood. Why do certain coordination patterns "
                "(in-phase/symmetric) come naturally while others (anti-phase) "
                "become unstable?"
            ),
            "task": (
                "Bimanual reaching movements with different amplitude and timing "
                "requirements. Behavioral analysis of temporal coupling between "
                "hands."
            ),
            "results": (
                "When the two hands had to move different distances, they tended to "
                "synchronize their movements temporally, with the shorter movement "
                "being slowed down. This seminal finding established that the motor "
                "system has a strong tendency toward temporal coupling between "
                "hands, consistent with the intrinsic coordination constraints "
                "that the manuscript's intrinsic neural correlations now provide "
                "a neural basis for."
            ),
        },
        {
            "title": "Perceptual basis of bimanual coordination",
            "authors": "Mechsner F, Kerzel D, Knoblich G, Prinz W",
            "year": "2001",
            "journal": "Nature, 414(6859), 69-73",
            "doi_or_link": "https://doi.org/10.1038/35102060",
            "relation_score": 7,
            "gap": (
                "Whether bimanual coordination is governed by perceptual/spatial "
                "(extrinsic) symmetry or by muscle-level (intrinsic) symmetry was "
                "debated. The dominant view favored motor/intrinsic constraints, "
                "but this study challenged that view."
            ),
            "task": (
                "Bimanual circle-drawing with visual transformations that decoupled "
                "spatial from muscular symmetry. Behavioral coordination stability."
            ),
            "results": (
                "Argued that bimanual coordination is governed primarily by "
                "perceptual (extrinsic) symmetry, not muscle homology. When visual "
                "feedback was transformed so that non-homologous muscle patterns "
                "produced spatially symmetric output, coordination was stable. This "
                "directly contrasts with the manuscript's finding of intrinsic "
                "(not extrinsic) neural organization, suggesting the brain's "
                "internal coding may differ from what behavioral performance implies."
            ),
        },
        {
            "title": "Muscle and movement representations in the primary motor cortex",
            "authors": "Kakei S, Hoffman DS, Strick PL",
            "year": "1999",
            "journal": "Science, 285(5436), 2136-2139",
            "doi_or_link": "https://doi.org/10.1126/science.285.5436.2136",
            "relation_score": 9,
            "gap": (
                "Whether individual M1 neurons encode movement direction in "
                "extrinsic (spatial) or intrinsic (joint/muscle) coordinates was "
                "debated. Different wrist postures can dissociate these frames, "
                "revealing the true coordinate system."
            ),
            "task": (
                "Monkey single-unit electrophysiology in M1 during wrist movements "
                "in different forearm postures (pronated vs. supinated) to "
                "dissociate extrinsic and intrinsic coordinates."
            ),
            "results": (
                "M1 neurons fall into distinct groups: some encode in extrinsic "
                "coordinates (direction tuning unchanged with posture), while "
                "others encode in intrinsic/muscle coordinates (tuning shifts with "
                "posture to follow muscle activation). This heterogeneity in M1 is "
                "directly relevant to the manuscript's question of whether "
                "bilateral representations are intrinsically or extrinsically "
                "organized. The wrist task is similar to the manuscript's."
            ),
        },
        {
            "title": "Primary motor cortex is involved in bimanual coordination",
            "authors": "Donchin O, Gribova A, Steinberg O, Bergman H, Vaadia E",
            "year": "1998",
            "journal": "Nature, 395(6699), 274-278",
            "doi_or_link": "https://doi.org/10.1038/26220",
            "relation_score": 8,
            "gap": (
                "During bimanual movements, it was unknown whether M1 neurons that "
                "respond to both hands do so in intrinsic (muscle) or extrinsic "
                "(spatial) coordinates. This is crucial for understanding "
                "bimanual interference at the neural level."
            ),
            "task": (
                "Monkey single-unit electrophysiology in M1 during bimanual "
                "reaching. Compared neural tuning for ipsilateral movements in "
                "intrinsic vs. extrinsic coordinates."
            ),
            "results": (
                "Many M1 neurons showed significant modulation during ipsilateral "
                "movements. This landmark study established that M1 neurons carry "
                "bilateral information and are actively involved in bimanual "
                "coordination, not just contralateral hand control. Laid the "
                "groundwork for testing whether ipsilateral representations are "
                "organized in intrinsic or extrinsic coordinates at the population "
                "level, as the manuscript does with fMRI."
            ),
        },
        {
            "title": "Integration of target and effector information in the human brain during reach planning",
            "authors": "Beurze SM, de Lange FP, Toni I, Medendorp WP",
            "year": "2007",
            "journal": "Journal of Neurophysiology, 97(1), 188-199",
            "doi_or_link": "https://doi.org/10.1152/jn.00456.2006",
            "relation_score": 5,
            "gap": (
                "How parietal cortex integrates target location (extrinsic) with "
                "effector selection (left/right hand) during reaching was unclear. "
                "Do parietal areas use extrinsic, intrinsic, or intermediate "
                "reference frames?"
            ),
            "task": (
                "Reaching with left or right hand to visual targets during fMRI. "
                "Manipulated effector and target location to dissociate reference "
                "frames in posterior parietal cortex."
            ),
            "results": (
                "Posterior parietal cortex showed integration of target and effector "
                "information during reach planning. The representations showed a "
                "gradient from extrinsic to intrinsic coordinates across the "
                "parietal-to-premotor pathway. This is consistent with the "
                "manuscript's finding that parietal areas show bilateral "
                "representations and may explain why SPL areas show the intrinsic "
                "organization found."
            ),
        },
        {
            "title": "Static and phasic cross-talk effects in discrete bimanual reversal movements",
            "authors": "Heuer H, Kleinsorge T, Spijkers W, Steglich C",
            "year": "2001",
            "journal": "Journal of Motor Behavior, 33(1), 67-85",
            "doi_or_link": "https://doi.org/10.1080/00222890109601904",
            "relation_score": 7,
            "gap": (
                "Previous work treated intrinsic and extrinsic constraints as "
                "competing accounts, but their cooperative interaction was "
                "unexplored. Do these two reference frames combine additively "
                "or interact during bimanual movements?"
            ),
            "task": (
                "Bimanual reversal movements manipulating both intrinsic and "
                "extrinsic congruency simultaneously. Behavioral analysis of "
                "cross-talk and interference patterns."
            ),
            "results": (
                "Intrinsic and extrinsic constraints cooperate but with "
                "under-additivity, meaning the combined effect is less than "
                "the sum of individual effects. The intrinsic constraint was "
                "dominant for initial movement phases. This under-additive "
                "combination parallels the manuscript's sub-additive neural "
                "findings during bimanual movements, where ipsilateral "
                "representations diminish."
            ),
        },

        # ===== SECTION 3: Ipsilateral Motor Representations =====
        {
            "title": "Ipsilateral motor cortex activity during unimanual hand movements relates to task complexity",
            "authors": "Verstynen T, Diedrichsen J, Albert N, Aparicio P, Ivry RB",
            "year": "2005",
            "journal": "Journal of Neurophysiology, 93(3), 1209-1222",
            "doi_or_link": "https://doi.org/10.1152/jn.00720.2004",
            "relation_score": 8,
            "gap": (
                "The functional significance of ipsilateral M1 activity during "
                "unimanual movements was debated. Was it a genuine motor "
                "representation, an artifact of incomplete movement isolation, "
                "or a mirror of the contralateral pattern?"
            ),
            "task": (
                "Unimanual finger sequences of varying complexity during fMRI. "
                "Compared ipsilateral activation as a function of task complexity."
            ),
            "results": (
                "Confirmed genuine ipsilateral representations in M1 that scale "
                "with task complexity. Ipsilateral activation increased with more "
                "complex sequences. This validates the manuscript's finding of "
                "significant ipsilateral encoding in M1 and supports that it "
                "is not an artifact but reflects real motor information that "
                "scales with motor demands."
            ),
        },
        {
            "title": "The cortical physiology of ipsilateral limb movements",
            "authors": "Bundy DT, Leuthardt EC",
            "year": "2019",
            "journal": "Trends in Neurosciences, 42(11), 825-839",
            "doi_or_link": "https://doi.org/10.1016/j.tins.2019.08.008",
            "relation_score": 7,
            "gap": (
                "A comprehensive framework for understanding ipsilateral motor "
                "cortex activity was lacking. What functional roles does "
                "ipsilateral M1 serve, and how should we interpret ipsilateral "
                "activation in health and disease?"
            ),
            "task": (
                "Review article synthesizing fMRI, ECoG, TMS, and patient studies "
                "on ipsilateral motor cortex contributions to limb control."
            ),
            "results": (
                "Ipsilateral M1 activity reflects multiple sources: postural "
                "stabilization, bilateral coordination, and potential backup "
                "pathways. The review provides a framework for interpreting "
                "the manuscript's finding that ipsilateral encoding is present "
                "during unimanual movements but reduced during bimanual "
                "movements, suggesting context-dependent gating of ipsilateral "
                "representations."
            ),
        },
        {
            "title": "The role of multiple contralesional motor areas for complex hand movements after internal capsular lesion",
            "authors": "Lotze M, Markert J, Sauseng P, Hoppe J, Plewnia C, Gerloff C",
            "year": "2006",
            "journal": "Journal of Neuroscience, 26(22), 6096-6102",
            "doi_or_link": "https://doi.org/10.1523/JNEUROSCI.4564-05.2006",
            "relation_score": 4,
            "gap": (
                "After stroke damaging contralateral motor pathways, patients "
                "sometimes recover using ipsilateral pathways. Whether "
                "ipsilateral motor areas genuinely support movement recovery was "
                "clinically important and speaks to the latent bilateral "
                "representations found in healthy brains."
            ),
            "task": (
                "Stroke patients with internal capsular lesions performed hand "
                "movements during fMRI. Compared ipsilateral motor area "
                "activation with healthy controls. Patient study."
            ),
            "results": (
                "Stroke patients showed dramatically increased contralesional "
                "(ipsilateral) motor cortex activation during movements of the "
                "affected hand. Multiple motor areas contributed, including M1, "
                "PMd, and SMA. Demonstrates that the ipsilateral representations "
                "documented in the manuscript (present but suppressed during "
                "bimanual movement) may serve as a latent backup system."
            ),
        },
        {
            "title": "Effector-independent motor sequence representations exist in extrinsic and intrinsic reference frames",
            "authors": "Wiestler T, Waters-Metenier S, Diedrichsen J",
            "year": "2014",
            "journal": "Journal of Neuroscience, 34(14), 5054-5064",
            "doi_or_link": "https://doi.org/10.1523/JNEUROSCI.5363-13.2014",
            "relation_score": 8,
            "gap": (
                "Whether ipsilateral movement representations in motor cortex "
                "are organized in intrinsic or extrinsic coordinates was unknown. "
                "If ipsilateral patterns share coordinate frames with contralateral "
                "patterns, they may reflect shared encoding principles."
            ),
            "task": (
                "Finger sequence learning with transfer tests across hands, measured "
                "with fMRI. Multivariate analysis of ipsilateral and contralateral "
                "sequence representations."
            ),
            "results": (
                "Motor sequence representations exist in both extrinsic and "
                "intrinsic reference frames, and these persist across effectors. "
                "Learning a sequence with one hand created representations "
                "accessible by the other hand, organized in both coordinate "
                "systems. This dual coding is relevant to interpreting the "
                "manuscript's finding of intrinsic ipsilateral encoding."
            ),
        },
        {
            "title": "Functional magnetic resonance imaging of motor cortex: hemispheric asymmetry and handedness",
            "authors": "Kim SG, Ashe J, Hendrich K, Ellermann JM, Merkle H, Ugurbil K, Georgopoulos AP",
            "year": "1993",
            "journal": "Science, 261(5121), 615-617",
            "doi_or_link": "https://doi.org/10.1126/science.8342027",
            "relation_score": 6,
            "gap": (
                "Whether the degree of hemispheric lateralization in M1 varies "
                "with handedness and whether ipsilateral M1 activation is present "
                "during simple hand movements was unknown in the early fMRI era."
            ),
            "task": (
                "Sequential finger-thumb opposition with dominant and non-dominant "
                "hands during fMRI. Compared contralateral vs. ipsilateral M1."
            ),
            "results": (
                "Found significant contralateral M1 activation for both hands, "
                "with additional ipsilateral M1 activation more prominent for the "
                "non-dominant hand. The ratio of contralateral to ipsilateral was "
                "about 3:1, remarkably similar to the 3.48x ratio reported in the "
                "manuscript. Established the baseline lateralization pattern that "
                "the manuscript builds upon."
            ),
        },

        # ===== SECTION 4: Bimanual Coordination - Monkey Electrophysiology =====
        {
            "title": "Neural interactions between motor cortical hemispheres during bimanual and unimanual arm movements",
            "authors": "Cardoso de Oliveira S, Gribova A, Donchin O, Bergman H, Vaadia E",
            "year": "2001",
            "journal": "European Journal of Neuroscience, 14(11), 1881-1896",
            "doi_or_link": "https://doi.org/10.1046/j.0953-816x.2001.01801.x",
            "relation_score": 8,
            "gap": (
                "The neural mechanism of bimanual interference in M1 was unknown. "
                "Do interhemispheric interactions between M1 populations change "
                "when the ipsilateral hand simultaneously moves, and could this "
                "cause interference?"
            ),
            "task": (
                "Monkey single-unit recording in bilateral M1 during unimanual and "
                "bimanual reaching. Analyzed interhemispheric neural correlations."
            ),
            "results": (
                "Interhemispheric correlations between M1 neurons changed "
                "significantly during bimanual compared to unimanual movements. "
                "This context-dependent modulation of inter-hemispheric "
                "communication could cause interference because the same "
                "population is influenced by the other hemisphere's activity. "
                "This parallels the manuscript's finding that ipsilateral encoding "
                "changes during bimanual movements."
            ),
        },
        {
            "title": "Neuronal activity in cortical motor areas related to ipsilateral, contralateral, and bilateral digit movements of the monkey",
            "authors": "Tanji J, Okano K, Sato KC",
            "year": "1988",
            "journal": "Journal of Neurophysiology, 60(1), 325-343",
            "doi_or_link": "https://doi.org/10.1152/jn.1988.60.1.325",
            "relation_score": 7,
            "gap": (
                "Whether single neurons in SMA and M1 are preferentially "
                "activated during bimanual vs. unimanual movements was unknown. "
                "SMA's specific role in bimanual coordination needed "
                "electrophysiological evidence."
            ),
            "task": (
                "Monkey single-unit recording in SMA and M1 during ipsilateral, "
                "contralateral, and bilateral digit movements."
            ),
            "results": (
                "Found neurons in SMA specifically active during bimanual "
                "movements and not during either unimanual movement alone, "
                "providing direct evidence for bimanual-specific representations. "
                "M1 neurons were primarily contralateral-hand driven with less "
                "bimanual specificity. This SMA specialization is consistent "
                "with the manuscript's finding that the interaction component "
                "was present across regions including SMA."
            ),
        },
        {
            "title": "Single-unit activity related to bimanual arm movements in the primary and supplementary motor cortices",
            "authors": "Donchin O, Gribova A, Steinberg O, Mitz AR, Bergman H, Vaadia E",
            "year": "2002",
            "journal": "Journal of Neurophysiology, 88(6), 3498-3517",
            "doi_or_link": "https://doi.org/10.1152/jn.00335.2001",
            "relation_score": 8,
            "gap": (
                "Whether M1 neurons carry bimanual-specific information beyond "
                "what can be predicted from their unimanual responses was "
                "unknown. If bimanual activity is simply the sum of unimanual "
                "components, M1 may not actively coordinate bimanual movements."
            ),
            "task": (
                "Monkey single-unit recording in M1 and SMA during unimanual and "
                "bimanual arm reaching. Compared actual bimanual responses to "
                "predictions from unimanual data."
            ),
            "results": (
                "Approximately 50% of M1 neurons showed bimanual responses "
                "that could not be predicted from their unimanual activities, "
                "demonstrating a genuine bimanual interaction at the single-"
                "neuron level. SMA showed even stronger bimanual-specific "
                "responses. This provides single-neuron evidence for the "
                "interaction component identified in the manuscript at the "
                "population level using fMRI."
            ),
        },
        {
            "title": "Neural activity in primary motor and dorsal premotor cortex in reaching tasks with the contralateral versus ipsilateral arm",
            "authors": "Cisek P, Crammond DJ, Kalaska JF",
            "year": "2003",
            "journal": "Journal of Neurophysiology, 89(2), 922-942",
            "doi_or_link": "https://doi.org/10.1152/jn.00607.2002",
            "relation_score": 7,
            "gap": (
                "The degree to which M1 and PMd neurons encode ipsilateral arm "
                "movements in primates was not systematically quantified. How "
                "lateralized are these areas, and does lateralization differ "
                "between M1 and premotor cortex?"
            ),
            "task": (
                "Monkey single-unit recording in M1 and PMd during unimanual "
                "reaching with contralateral vs. ipsilateral arm."
            ),
            "results": (
                "Both M1 and PMd showed significant modulation during "
                "ipsilateral arm movements, though weaker than contralateral. "
                "PMd showed more bilateral representation than M1, with "
                "similar directional tuning for both arms. The ipsilateral "
                "representations in PMd were consistent with intrinsic "
                "coordinates. This M1-PMd lateralization gradient matches "
                "the manuscript's fMRI finding."
            ),
        },
        {
            "title": "Neural population dynamics during reaching",
            "authors": "Churchland MM, Cunningham JP, Kaufman MT, Foster JD, Nuyujukian P, Ryu SI, Shenoy KV",
            "year": "2012",
            "journal": "Nature, 487(7405), 51-56",
            "doi_or_link": "https://doi.org/10.1038/nature11129",
            "relation_score": 3,
            "gap": (
                "Motor cortex activity during reaching was traditionally analyzed "
                "with static tuning curves, but a dynamics-based view was needed. "
                "How do M1 neural populations evolve over time during reaching?"
            ),
            "task": (
                "Monkey single-unit recording in M1 and PMd during delayed "
                "center-out reaching. Population dynamics analyzed with jPCA "
                "dimensionality reduction."
            ),
            "results": (
                "M1 population activity during reaching follows rotational "
                "dynamics consistent across conditions. This dynamical systems "
                "view provides context for understanding how the same M1 "
                "population might handle bilateral representations: the dynamics "
                "may occupy orthogonal subspaces for each hand, consistent "
                "with the manuscript's decomposition approach."
            ),
        },
        {
            "title": "Motor cortex signals for each arm are mixed across hemispheres and neurons yet partitioned within the population response",
            "authors": "Ames KC, Churchland MM",
            "year": "2019",
            "journal": "eLife, 8, e46159",
            "doi_or_link": "https://doi.org/10.7554/eLife.46159",
            "relation_score": 9,
            "gap": (
                "If M1 represents both hands, how does it avoid cross-talk? "
                "One hypothesis is that left and right hand information occupies "
                "orthogonal neural subspaces, preventing interference between "
                "the two hand representations."
            ),
            "task": (
                "Monkey electrophysiology in M1 during unimanual and bimanual "
                "reaching. Dimensionality reduction to identify hand-specific "
                "subspaces."
            ),
            "results": (
                "Individual neurons were mixed in their selectivity for each arm, "
                "but at the population level, left and right arm representations "
                "occupied largely orthogonal neural subspaces. During bimanual "
                "movements, each hand's representation remained in its "
                "respective subspace, minimizing cross-talk. This orthogonality "
                "mechanism is directly relevant to the manuscript's discussion "
                "of how contralateral encoding remains stable while ipsilateral "
                "encoding disappears during bimanual movements."
            ),
        },
        {
            "title": "Reorganization between preparatory and movement population responses in motor cortex",
            "authors": "Elsayed GF, Lara AH, Kaufman MT, Churchland MM, Cunningham JP",
            "year": "2016",
            "journal": "Nature Communications, 7, 13239",
            "doi_or_link": "https://doi.org/10.1038/ncomms13239",
            "relation_score": 4,
            "gap": (
                "How motor cortex separates movement preparation from execution "
                "within the same neural population was unclear. The null-space "
                "hypothesis proposed that preparatory activity resides in a "
                "subspace that does not drive motor output."
            ),
            "task": (
                "Monkey electrophysiology in M1 and PMd during delayed reaching. "
                "Subspace analysis separating preparatory from movement activity."
            ),
            "results": (
                "Preparatory and movement-related activity occupied orthogonal "
                "neural subspaces within the same population. This orthogonal "
                "subspace mechanism is the same principle that may allow M1 to "
                "represent both hands without interference. Relevant to "
                "understanding how the ipsilateral representation can exist "
                "without causing involuntary movements of the non-moving hand."
            ),
        },

        # ===== SECTION 5: SMA, Premotor, & Parietal Cortex =====
        {
            "title": "Supplementary motor area of the monkey's cerebral cortex: short- and long-term deficits after unilateral ablation and the effects of subsequent callosal section",
            "authors": "Brinkman C",
            "year": "1984",
            "journal": "Journal of Neuroscience, 4(4), 918-929",
            "doi_or_link": "https://doi.org/10.1523/JNEUROSCI.04-04-00918.1984",
            "relation_score": 5,
            "gap": (
                "Whether SMA damage specifically impairs bimanual movements "
                "while sparing unimanual ones was unknown. Causal evidence for "
                "SMA's unique bimanual role was needed."
            ),
            "task": (
                "Monkey lesion study: unilateral SMA ablation followed by testing "
                "unimanual and bimanual hand movements, then callosal section. "
                "Behavioral assessment."
            ),
            "results": (
                "SMA lesions caused severe deficits in bimanual coordination "
                "while unimanual movements remained relatively intact. Animals "
                "could still move each hand independently but could not "
                "coordinate them together. Subsequent callosal section worsened "
                "bimanual deficits. This provides causal evidence that SMA is "
                "necessary for bimanual coordination, supporting the manuscript's "
                "finding of bimanual-specific activity in SMA."
            ),
        },
        {
            "title": "Clinical consequences of corticectomies involving the supplementary motor area in man",
            "authors": "Laplane D, Talairach J, Meininger V, Bancaud J, Orgogozo JM",
            "year": "1977",
            "journal": "Journal of the Neurological Sciences, 34(3), 301-314",
            "doi_or_link": "https://doi.org/10.1016/0022-510x(77)90148-4",
            "relation_score": 5,
            "gap": (
                "Whether SMA damage in humans causes bimanual-specific deficits "
                "similar to monkey lesion studies was needed for cross-species "
                "validation of SMA's role."
            ),
            "task": (
                "Neurological examination of patients with corticectomies involving "
                "the SMA. Assessment of unimanual vs. bimanual motor tasks. "
                "Patient study."
            ),
            "results": (
                "Patients with SMA lesions showed bimanual incoordination: "
                "normal unimanual movements but difficulty with bimanual tasks, "
                "especially asymmetric (non-mirror) movements. This human "
                "patient evidence complements monkey studies and is consistent "
                "with the manuscript's finding that SMA carries bimanual-specific "
                "interaction information."
            ),
        },
        {
            "title": "The role of anterior cingulate cortex and precuneus in the coordination of motor behaviour",
            "authors": "Wenderoth N, Debaere F, Sunaert S, Swinnen SP",
            "year": "2005",
            "journal": "European Journal of Neuroscience, 22(1), 235-246",
            "doi_or_link": "https://doi.org/10.1111/j.1460-9568.2005.04176.x",
            "relation_score": 7,
            "gap": (
                "Which cortical network supports bimanual coordination, and "
                "how does network activity change with coordination complexity? "
                "The specific roles of cingulate cortex and precuneus were unknown."
            ),
            "task": (
                "Bimanual wrist movements with different coordination patterns "
                "(in-phase, anti-phase, multi-frequency ratios) during fMRI. "
                "Univariate activation mapping."
            ),
            "results": (
                "Bimanual coordination recruited anterior cingulate cortex, SMA, "
                "PMd, and precuneus beyond unimanual activation. More complex "
                "coordination patterns increased activation in these regions "
                "specifically. M1 activation was primarily contralateral. "
                "This network mapping is consistent with the manuscript's ROI "
                "selection and finding that premotor and medial wall areas show "
                "bimanual-specific activity."
            ),
        },
        {
            "title": "Cerebellar and premotor function in bimanual coordination: parametric neural responses to spatiotemporal complexity and cycling frequency",
            "authors": "Debaere F, Wenderoth N, Sunaert S, van Hecke P, Swinnen SP",
            "year": "2004",
            "journal": "NeuroImage, 21(4), 1416-1427",
            "doi_or_link": "https://doi.org/10.1016/j.neuroimage.2003.12.011",
            "relation_score": 7,
            "gap": (
                "The neural substrates distinguishing stable (in-phase) from "
                "unstable (anti-phase) bimanual coordination were not well "
                "characterized, particularly the cerebellar contribution."
            ),
            "task": (
                "Bimanual cyclical wrist/forearm movements in in-phase, "
                "anti-phase, and 90-degree out-of-phase patterns during fMRI."
            ),
            "results": (
                "Anti-phase and complex coordination patterns produced greater "
                "activation in SMA, cerebellum, PMd, and secondary somatosensory "
                "cortex compared to in-phase. M1 showed similar activation "
                "for all patterns. This suggests the additional neural cost "
                "of non-intrinsic coordination is borne by premotor, "
                "supplementary, and cerebellar areas, consistent with the "
                "manuscript's interaction component findings."
            ),
        },
        {
            "title": "Distribution of eye- and arm-movement-related neuronal activity in the SEF and in the SMA and pre-SMA of monkeys",
            "authors": "Fujii N, Mushiake H, Tanji J",
            "year": "2002",
            "journal": "Journal of Neurophysiology, 87(4), 2158-2166",
            "doi_or_link": "https://doi.org/10.1152/jn.00867.2001",
            "relation_score": 3,
            "gap": (
                "How SMA and pre-SMA neurons are organized with respect to "
                "arm and eye movement control was unknown. Functional "
                "segregation within these areas could explain their roles in "
                "coordinating multi-effector and bimanual movements."
            ),
            "task": (
                "Monkey single-unit recording in SMA, pre-SMA, and SEF during "
                "arm and eye movements. Analyzed spatial distribution of "
                "movement-related neurons."
            ),
            "results": (
                "SMA contained a high proportion of arm-movement-related "
                "neurons organized somatotopically, while pre-SMA contained "
                "more mixed arm/eye neurons. This functional organization "
                "supports SMA's role in limb coordination, including bimanual "
                "control, consistent with the manuscript's finding of "
                "bimanual-specific activity in medial motor areas."
            ),
        },
        {
            "title": "Coordination of uncoupled bimanual movements by strictly timed interhemispheric connectivity",
            "authors": "Liuzzi G, Horniss V, Zimerman M, Gerloff C, Hummel FC",
            "year": "2011",
            "journal": "Journal of Neuroscience, 31(25), 9111-9117",
            "doi_or_link": "https://doi.org/10.1523/JNEUROSCI.0046-11.2011",
            "relation_score": 8,
            "gap": (
                "The temporal dynamics of interhemispheric communication during "
                "bimanual coordination were unknown. Whether precisely timed "
                "connectivity between motor cortices is required for coordination "
                "needed causal testing."
            ),
            "task": (
                "Paired-pulse TMS between motor cortices during bimanual finger "
                "movements. Measured temporal specificity of interhemispheric "
                "inhibition/facilitation."
            ),
            "results": (
                "Precisely timed interhemispheric connectivity was required for "
                "successful bimanual coordination. Disrupting this timing with "
                "TMS impaired coordination without affecting unimanual movements. "
                "This causal evidence supports the idea that the intrinsic "
                "coupling between hemispheric representations in the manuscript "
                "depends on callosal communication with precise timing."
            ),
        },

        # ===== SECTION 6: Behavioral Foundations & Context =====
        {
            "title": "Bimanual movement control: information processing and interaction effects",
            "authors": "Marteniuk RG, MacKenzie CL, Baba DM",
            "year": "1984",
            "journal": "Quarterly Journal of Experimental Psychology, 36A(2), 335-365",
            "doi_or_link": "https://doi.org/10.1080/14640748408402163",
            "relation_score": 6,
            "gap": (
                "Whether bimanual interference (cross-talk) operates at the "
                "level of spatial plans or motor commands was unknown. "
                "Dissociating spatial from motoric interference reveals the "
                "locus of interlimb coupling."
            ),
            "task": (
                "Bimanual aiming movements to targets of different sizes "
                "and distances. Behavioral analysis of movement time and "
                "accuracy for each hand."
            ),
            "results": (
                "When the two hands moved to targets of different difficulty, "
                "each hand's movement time was pulled toward the other. This "
                "assimilation effect demonstrated strong interlimb coupling at "
                "the motor planning level, consistent with shared "
                "representations. The behavioral interference matches the "
                "slight bimanual interference in the manuscript for "
                "incongruent conditions."
            ),
        },
        {
            "title": "Callosotomy patients exhibit temporal uncoupling during continuous bimanual movements",
            "authors": "Kennerley SW, Diedrichsen J, Hazeltine E, Semjen A, Ivry RB",
            "year": "2002",
            "journal": "Nature Neuroscience, 5(4), 376-381",
            "doi_or_link": "https://doi.org/10.1038/nn822",
            "relation_score": 8,
            "gap": (
                "Whether the corpus callosum is necessary for temporal "
                "coupling between hands during bimanual movements was debated. "
                "Some temporal coupling might be mediated subcortically."
            ),
            "task": (
                "Callosotomy patients and controls performing continuous "
                "bimanual rhythmic movements. Measured temporal coupling. "
                "Patient study."
            ),
            "results": (
                "Callosotomy patients showed temporal uncoupling during "
                "continuous (but not discrete) bimanual movements, demonstrating "
                "that the corpus callosum is specifically necessary for ongoing "
                "temporal coordination. Discrete movements retained some "
                "subcortical coupling. The intrinsic coupling found in the "
                "manuscript at the cortical level likely depends on this "
                "callosal transmission."
            ),
        },
        {
            "title": "Intermanual coordination: from behavioural principles to neural-network interactions",
            "authors": "Swinnen SP",
            "year": "2002",
            "journal": "Nature Reviews Neuroscience, 3(5), 348-359",
            "doi_or_link": "https://doi.org/10.1038/nrn807",
            "relation_score": 7,
            "gap": (
                "A comprehensive review was needed to integrate behavioral "
                "principles of bimanual coordination with emerging "
                "neuroimaging and electrophysiology findings."
            ),
            "task": (
                "Review covering behavioral experiments, fMRI, PET, lesion "
                "studies, and electrophysiology on bimanual coordination."
            ),
            "results": (
                "Identified key principles: bimanual coordination is "
                "constrained by both spatial and muscular symmetry; SMA and "
                "premotor areas are critical; the corpus callosum enables "
                "interhemispheric communication for asymmetric movements; M1 "
                "is primarily contralateral but carries ipsilateral "
                "information. Provides the conceptual framework within which "
                "the manuscript's findings should be interpreted."
            ),
        },
        {
            "title": "Two hands, one brain: cognitive neuroscience of bimanual skill",
            "authors": "Swinnen SP, Wenderoth N",
            "year": "2004",
            "journal": "Trends in Cognitive Sciences, 8(1), 18-25",
            "doi_or_link": "https://doi.org/10.1016/j.tics.2003.10.017",
            "relation_score": 6,
            "gap": (
                "How intermanual interactions during planning vs. execution "
                "contribute to bimanual interference needed synthesis across "
                "behavioral and imaging modalities."
            ),
            "task": (
                "Review of behavioral and neuroimaging studies on bimanual "
                "skill learning and interference mechanisms."
            ),
            "results": (
                "Interference arises during both planning and execution, "
                "but the type differs. Planning-phase interference reflects "
                "shared representations in premotor/parietal areas (consistent "
                "with the manuscript's bilateral premotor encoding), while "
                "execution-phase interference involves M1 and subcortical "
                "structures. The intrinsic constraint is more dominant during "
                "execution, matching the intrinsic organization found."
            ),
        },
        {
            "title": "Pattern component modeling: a flexible approach for understanding the representational structure of brain activity patterns",
            "authors": "Diedrichsen J, Yokoi A, Arbuckle SA",
            "year": "2018",
            "journal": "NeuroImage, 180(Pt A), 119-133",
            "doi_or_link": "https://doi.org/10.1016/j.neuroimage.2017.08.051",
            "relation_score": 6,
            "gap": (
                "How to appropriately model and test the representational "
                "structure of brain activity patterns needed a flexible "
                "statistical framework beyond RSA."
            ),
            "task": (
                "Methodological paper with finger movement fMRI examples. "
                "Introduces pattern component modeling (PCM) for testing "
                "representational models."
            ),
            "results": (
                "Introduced PCM as a likelihood-based framework for testing "
                "hypotheses about representational structure, complementing "
                "the crossnobis dissimilarity used in the manuscript. "
                "Cross-validated measures are essential to avoid bias, and "
                "component models can test specific hypotheses about how "
                "bimanual patterns decompose into additive vs. interactive "
                "components."
            ),
        },
        {
            "title": "On the relations between the direction of two-dimensional arm movements and cell discharge in primate motor cortex",
            "authors": "Georgopoulos AP, Kalaska JF, Caminiti R, Massey JT",
            "year": "1982",
            "journal": "Journal of Neuroscience, 2(11), 1527-1537",
            "doi_or_link": "https://doi.org/10.1523/JNEUROSCI.02-11-01527.1982",
            "relation_score": 3,
            "gap": (
                "How M1 neurons encode movement direction during reaching was "
                "a foundational question. Establishing directional tuning was "
                "needed for all subsequent work on contralateral and "
                "ipsilateral direction encoding."
            ),
            "task": (
                "Monkey single-unit recording in M1 during center-out "
                "reaching to 8 targets. Characterized directional tuning "
                "curves and population vector decoding."
            ),
            "results": (
                "M1 neurons are broadly tuned for movement direction, each "
                "with a preferred direction. The population vector accurately "
                "predicted reach direction. This foundational work on "
                "directional encoding underpins the manuscript's analysis of "
                "directional encoding in motor regions. The same center-out "
                "paradigm is adapted in the manuscript for wrist movements."
            ),
        },
        {
            "title": "Structure of population activity in primary motor cortex for single finger flexion and extension",
            "authors": "Arbuckle SA, Weiler J, Kirk EA, Rice CL, Schieber M, Pruszynski JA, Ejaz N, Diedrichsen J",
            "year": "2020",
            "journal": "Journal of Neuroscience, 40(48), 9210-9223",
            "doi_or_link": "https://doi.org/10.1523/JNEUROSCI.0999-20.2020",
            "relation_score": 5,
            "gap": (
                "How the hand area of primary motor cortex (M1) is organized to control "
                "fine finger movements remains unclear, especially with respect to "
                "opposite movement directions. Because flexion and extension of the same "
                "finger never co-occur, a natural-statistics account would predict "
                "distinct representations—yet it was unknown whether population activity "
                "patterns in M1 separate flexion from extension, and how conclusions "
                "depend on measurement modality (fMRI vs spiking)."
            ),
            "task": (
                "Single-finger near-isometric flexion and extension task. Measured "
                "multivariate activity patterns with 7T fMRI in humans and compared them "
                "to neural spiking patterns from two monkeys performing the identical "
                "task; also related findings to muscular activity patterns."
            ),
            "results": (
                "Human fMRI patterns were distinct across different fingers, but were "
                "surprisingly similar for flexion versus extension of the same finger. "
                "In contrast, monkey M1 spiking patterns were distinct for both finger "
                "identity and movement direction (similar to muscular patterns). The "
                "fMRI–spiking discrepancy is consistent with an architecture where "
                "populations controlling flexion and extension of the same finger "
                "produce distinct outputs but share inputs and local recurrent "
                "processing and/or are spatially intermingled at a scale that makes "
                "opponent directions appear similar in coarse-scale fMRI. This provides "
                "testable hypotheses about how M1 organizes fine finger control."
            ),
        },
        {
            "title": "Learned dynamics of reaching movements generalize from dominant to nondominant arm",
            "authors": "Criscimagna-Hemminger SE, Donchin O, Gazzaniga MS, Shadmehr R",
            "year": "2003",
            "journal": "Journal of Neurophysiology, 89(1), 168-176",
            "doi_or_link": "https://doi.org/10.1152/jn.00622.2002",
            "relation_score": 4,
            "gap": (
                "Whether motor learning transfers between limbs and in which "
                "coordinate system was unknown. Transfer direction reveals the "
                "coordinate system of shared motor representations."
            ),
            "task": (
                "Force-field adaptation with the dominant arm, then tested with "
                "the nondominant arm. Behavioral analysis of interlimb transfer."
            ),
            "results": (
                "Learning transferred from dominant to nondominant arm, and the "
                "pattern of transfer was consistent with generalization in "
                "extrinsic coordinates. This provided key evidence about how "
                "internal models are shared across limbs. Relevant to the "
                "manuscript's question of how representations relate between "
                "hands, whether in intrinsic or extrinsic coordinates."
            ),
        },
        {
            "title": "Reach plans in eye-centered coordinates",
            "authors": "Batista AP, Buneo CA, Snyder LH, Andersen RA",
            "year": "1999",
            "journal": "Science, 285(5425), 257-260",
            "doi_or_link": "https://doi.org/10.1126/science.285.5425.257",
            "relation_score": 2,
            "gap": (
                "In what coordinate frame does parietal cortex represent "
                "reach targets? Whether parietal reach plans are coded in "
                "eye-centered (extrinsic) or hand/body-centered (intrinsic) "
                "coordinates was debated."
            ),
            "task": (
                "Monkey electrophysiology in parietal reach region during "
                "reaches with different eye and hand positions to dissociate "
                "reference frames."
            ),
            "results": (
                "Parietal reach region neurons encoded reach targets in "
                "eye-centered coordinates, an extrinsic frame of reference. "
                "This suggests that parietal cortex encodes goals in sensory "
                "coordinates before transformation into motor coordinates. "
                "Relevant to understanding the coordinate transformation "
                "gradient from extrinsic to intrinsic across the parietal-"
                "to-premotor pathway that the manuscript's intrinsic "
                "organization reflects."
            ),
        },
        {
            "title": "Motor cortical representation of speed and direction during reaching",
            "authors": "Moran DW, Schwartz AB",
            "year": "1999",
            "journal": "Journal of Neurophysiology, 82(5), 2676-2692",
            "doi_or_link": "https://doi.org/10.1152/jn.1999.82.5.2676",
            "relation_score": 2,
            "gap": (
                "Whether M1 neurons encode movement kinematics (direction, speed) "
                "or dynamics (force, torque) during reaching was debated. This "
                "distinction determines whether M1 operates in extrinsic or "
                "intrinsic coordinates."
            ),
            "task": (
                "Monkey electrophysiology in M1 during continuous drawing "
                "movements. Compared neural correlates with kinematic (speed, "
                "direction) and dynamic (force) variables."
            ),
            "results": (
                "M1 neurons encoded both speed and direction of hand movement. "
                "The population showed a strong relationship to movement "
                "kinematics. This mixed representation means M1 does not "
                "purely operate in one coordinate frame, which helps explain "
                "why the manuscript finds intrinsic organization at the "
                "population level even though single neurons show mixed coding."
            ),
        },
    ]


# ------------------------------------------------------------------ #
#  Main
# ------------------------------------------------------------------ #

def main():
    pdf = LitReviewPDF()
    pdf.alias_nb_pages()
    pdf.set_auto_page_break(auto=True, margin=20)
    pdf.add_title_page()

    papers = get_papers()

    sections = [
        ("Section 1: Bimanual Activity Decomposition & Neural Encoding (Papers 1-10)", 0, 10),
        ("Section 2: Intrinsic vs. Extrinsic Reference Frames (Papers 11-17)", 10, 17),
        ("Section 3: Ipsilateral Motor Representations (Papers 18-22)", 17, 22),
        ("Section 4: Bimanual Coordination - Monkey Electrophysiology (Papers 23-31)", 22, 31),
        ("Section 5: SMA, Premotor, & Parietal Cortex (Papers 32-38)", 31, 38),
        ("Section 6: Behavioral Foundations & Context (Papers 39-45)", 38, 45),
    ]

    for section_title, start, end in sections:
        pdf.add_page()
        pdf.add_section_header(section_title)
        for i in range(start, min(end, len(papers))):
            pdf.add_paper(i + 1, **papers[i])

    output_path = "/Users/aghavamp/Desktop/Projects/bimanual_wrist/literature_review_bimanual.pdf"
    pdf.output(output_path)
    print(f"PDF saved to: {output_path}")
    print(f"Total papers: {len(papers)}")


if __name__ == "__main__":
    main()
