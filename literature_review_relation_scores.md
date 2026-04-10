# Literature review (relation scores)

This file was generated from `tools/generate_lit_review.py`.

## Verification status

- **Directly verified against provided PDFs**: papers 1, 3, 5, 7, 8, 42.

- **Not fully PDF-verified in this pass**: the remaining papers were checked for bibliographic consistency (title/authors/year/journal/DOI) and topical plausibility, but the `gap/task/results` text was not line-by-line verified from the full paper PDF here.


## Papers

### 1. Neural organization of hierarchical motor sequence representations in the human neocortex

- **Authors**: Yokoi A, Diedrichsen J

- **Year**: 2019

- **Journal**: Neuron, 103(6), 1178-1190.e7

- **Link/DOI**: https://doi.org/10.1016/j.neuron.2019.06.017

- **Relation score (1–10)**: 2

- **Verification**: PDF-verified


**Gap**

Although it is widely accepted that movement sequences are represented hierarchically (individual movements, chunks, and whole sequences), it was unclear how these different levels map onto human neocortex. In particular, do different cortical regions represent distinct levels of the hierarchy (anatomical separation), or do chunk and sequence representations coexist within the same regions?


**Task & Modality**

Participants learned and executed 8 trained right-hand finger-press sequences explicitly organized into four chunks (2–3 presses each). Used representational fMRI analyses (multivariate pattern similarity) to test for representations of individual finger presses, chunks, and whole sequences across cortex.


**Key results**

Found clear evidence for all three representational levels. Anatomically, individual movements were uniquely represented in primary motor cortex (M1), whereas movement chunks and entire sequences were jointly represented with substantial overlap in premotor and parietal cortices. These findings challenge a strict anatomical hierarchy where each level is cleanly separated into different regions, and instead highlight a special distinction between representations of individual movements and sequential context.


---

### 2. Effector-invariant movement encoding in the human motor system

- **Authors**: Haar S, Dinstein I, Shelef I, Donchin O

- **Year**: 2017

- **Journal**: Journal of Neuroscience, 37(37), 9054-9063

- **Link/DOI**: https://doi.org/10.1523/JNEUROSCI.1663-17.2017

- **Relation score (1–10)**: 9

- **Verification**: Not PDF-verified


**Gap**

Whether ipsilateral motor regions encode movements in the same coordinate system as contralateral movements (extrinsic/spatial vs. intrinsic/joint-based) was unclear. Understanding this relationship reveals how the brain manages bilateral motor control without causing interference.


**Task & Modality**

Center-out reaching with each arm to 4 targets, measured with fMRI and multivariate pattern analysis (crossnobis dissimilarity).


**Key results**

Ipsilateral and contralateral movements with symmetric joint configurations were encoded similarly by neural populations in M1, PMd, SMA, and SPL, providing evidence for effector-invariant encoding in intrinsic/joint coordinates. This is directly relevant to the manuscript's finding that contra- and ipsilateral activity patterns are intrinsically correlated rather than extrinsically correlated. The study used the same crossnobis metric as the manuscript.


---

### 3. Future movement plans interact in sequential arm movements

- **Authors**: Kashefi M, Reschechtko S, Ariani G, Shahbazi M, Tan A, Diedrichsen J, Pruszynski JA

- **Year**: 2024

- **Journal**: eLife, 13, e94485

- **Link/DOI**: https://doi.org/10.7554/eLife.94485

- **Relation score (1–10)**: 2

- **Verification**: PDF-verified


**Gap**

Real-world actions require planning future movements while executing the current one. It was unclear (1) how many future reaches can be planned simultaneously in longer sequences and (2) whether planning processes for multiple future movements are independent or interact with each other and with ongoing control.


**Task & Modality**

Behavioral continuous sequential reaching: participants executed 14-reach sequences in a planar robotic exoskeleton while the number of visible future targets (planning horizon H1–H5) and dwell times (75/200/400 ms) were manipulated. Used target-jump perturbations to probe whether +2 reaches are planned and whether future plans interact.


**Key results**

Participants planned at least two future reaches during execution of the current reach. Jumping the +2 target revealed partial commitment to the pre-jump +2 location, indicating advance planning of +2. Critically, planning processes interacted: correcting an ongoing reach to a +1 target jump was slower/longer when participants could also plan +2 (H3) compared to when only +1 could be planned (H2). Additionally, current reach curvature was influenced by the next reach only when the planning processes temporally overlapped, consistent with resource sharing or interaction between future movement plans.


---

### 4. Where one hand meets the other: limb-specific and action-dependent movement plans decoded from preparatory signals in single human frontoparietal brain areas

- **Authors**: Gallivan JP, McLean DA, Flanagan JR, Culham JC

- **Year**: 2013

- **Journal**: Journal of Neuroscience, 33(5), 1991-2008

- **Link/DOI**: https://doi.org/10.1523/JNEUROSCI.0541-12.2013

- **Relation score (1–10)**: 5

- **Verification**: Not PDF-verified


**Gap**

While contralateral limb planning is well-established in frontoparietal cortex, the degree to which ipsilateral limb actions are represented was underappreciated. Can fMRI multivariate methods decode which hand will be used and what action will be performed from preparatory activity?


**Task & Modality**

Reaching vs. grasping with left or right hand while maintaining fixation. Decoded from fMRI preparatory signals using MVPA.


**Key results**

Found much stronger ipsilateral limb representation than expected. Subregions of posterior parietal cortex, PMd, PMv, pre-SMA, and even M1 could decode upcoming ipsilateral hand actions. This challenges the strict contralateral view and supports the manuscript's finding that ipsilateral encoding is present even in M1 and S1. Premotor and parietal areas showed approximately symmetric bilateral representations.


---

### 5. Skill learning strengthens cortical representations of motor sequences

- **Authors**: Wiestler T, Diedrichsen J

- **Year**: 2013

- **Journal**: eLife, 2, e00801

- **Link/DOI**: https://doi.org/10.7554/eLife.00801

- **Relation score (1–10)**: 3

- **Verification**: PDF-verified


**Gap**

Whether multivariate fMRI methods can reveal fine-grained movement representations in motor cortex, beyond what univariate activation shows, and whether learning strengthens these representations needed validation.


**Task & Modality**

Finger sequence learning with fMRI. Participants trained for 4 days on four left-hand finger sequences and then performed trained vs. untrained sequences during scanning. Used multivariate pattern analysis (MVPA) to discriminate sequence-specific activity patterns and quantify how reliably sequences could be distinguished.


**Key results**

Motor-skill acquisition was associated with the emergence of more distinguishable sequence-specific activity patterns for trained sequences, even without increases in spatially averaged activity. Both trained and untrained sequences could be discriminated in primary and secondary motor areas, but trained sequences were classified more reliably (especially in SMA). This demonstrates that learning can strengthen representational specificity even when univariate activation does not increase.


---

### 6. Distinct representation of ipsilateral hand movements in sensorimotor areas

- **Authors**: Bruurmijn MLCM, Raemaekers M, Branco MP, Ramsey NF, Vansteensel MJ

- **Year**: 2021

- **Journal**: European Journal of Neuroscience, 54(10), 7599-7608

- **Link/DOI**: https://doi.org/10.1111/ejn.15501

- **Relation score (1–10)**: 6

- **Verification**: Not PDF-verified


**Gap**

While contralateral motor cortex involvement is well-established, the role of ipsilateral sensorimotor cortex during unilateral movements remained debated. Can different ipsilateral hand gestures be decoded, and from which sub-regions of sensorimotor cortex?


**Task & Modality**

Multiple hand gestures during 7T fMRI. Multivariate decoding of ipsilateral vs. contralateral hand movement identity.


**Key results**

Ipsilateral hand gestures could be distinguished from each other using multivariate pattern analysis. Ipsilateral activation was strongest in anterior precentral gyrus and posterior postcentral gyrus. Critically, ipsilateral patterns were distinct from contralateral patterns, suggesting the same cortical region uses different population codes for each hand. Relevant to the manuscript's finding of ipsilateral encoding across motor regions.


---

### 7. The role of human primary motor cortex in the production of skilled finger sequences

- **Authors**: Yokoi A, Arbuckle SA, Diedrichsen J

- **Year**: 2018

- **Journal**: Journal of Neuroscience, 38(6), 1430-1442

- **Link/DOI**: https://doi.org/10.1523/JNEUROSCI.2798-17.2017

- **Relation score (1–10)**: 4

- **Verification**: PDF-verified


**Gap**

Although M1 is essential for producing individuated finger movements, it was unclear whether M1 also represents learned sequences of multiple finger movements (i.e., sequential context), or whether sequence representations reside primarily in premotor/parietal areas with M1 reflecting only execution of component finger presses.


**Task & Modality**

Participants practiced finger press sequences and performed them during 3T and 7T fMRI. Used multivariate representational analyses to test for sequence representations in M1 versus premotor/parietal cortex. Compared multi-finger sequence patterns to patterns for constituent single-finger movements, and used passive replay to test whether sequence-pattern effects reflect execution rather than hemodynamic nonlinearities.


**Key results**

After intensive practice, premotor and parietal areas encoded the different movement sequences, but there was little or no evidence for sequence representations in M1. Instead, M1 activity patterns during sequences were fully explained by a linear combination of patterns for individual finger movements, with a disproportionately strong weight on the first finger of the sequence (the 'first-finger effect'). Passive replay showed this effect is linked to active execution processes rather than fMRI hemodynamic nonlinearity. The results suggest M1 primarily reflects execution of component presses, while sequence context is represented in premotor/parietal cortex.


---

### 8. Two distinct ipsilateral cortical representations for individuated finger movements

- **Authors**: Diedrichsen J, Wiestler T, Krakauer JW

- **Year**: 2013

- **Journal**: Cerebral Cortex, 23(6), 1362-1377

- **Link/DOI**: https://doi.org/10.1093/cercor/bhs120

- **Relation score (1–10)**: 10

- **Verification**: PDF-verified


**Gap**

Ipsilateral motor cortex shows activity changes during unimanual movements, but its functional significance was unclear. Specifically, do ipsilateral cortical areas carry finger-specific information, and does the nature of ipsilateral representations change between unimanual and bimanual contexts (i.e., when both hemispheres are actively engaged)?


**Task & Modality**

Human fMRI with multivoxel pattern analysis of fine-grained activation patterns during paced isometric finger presses. Experiment 1 tested unimanual left vs. right finger presses (all 10 digits). Experiment 2 tested unimanual and bimanual presses (digits 1, 3, 5) to examine how contra- and ipsilateral finger representations interact during bimanual actions.


**Key results**

Showed two fundamentally different ipsilateral representations. During unimanual ipsilateral presses, primary sensory and motor cortices contained finger-specific patterns that were nearly identical to contralateral mirror-symmetric patterns, despite overall suppression; this mirrored component vanished during bimanual actions. A second ipsilateral representation emerged during bimanual actions in caudal premotor and anterior parietal cortices, where ipsilateral actions were encoded as a nonlinear modulation of patterns related to contralateral actions. Together, the results show that ipsilateral representations change their informational content and likely functional role depending on behavioral context.


---

### 9. Hand use predicts the structure of representations in sensorimotor cortex

- **Authors**: Ejaz N, Hamada M, Diedrichsen J

- **Year**: 2015

- **Journal**: Nature Neuroscience, 18(7), 1034-1040

- **Link/DOI**: https://doi.org/10.1038/nn.4038

- **Relation score (1–10)**: 5

- **Verification**: Not PDF-verified


**Gap**

What determines the structure of movement representations in sensorimotor cortex? Is it anatomy, biomechanics, or the statistics of natural hand use? Understanding this is essential for interpreting the intrinsic organization found between hands.


**Task & Modality**

Single-finger presses measured with fMRI. RSA comparing brain patterns to models based on muscle anatomy and natural hand-use statistics.


**Key results**

The structure of finger representations in M1 was best predicted by the statistics of natural hand use (how often fingers are used together in daily life), not by muscle anatomy alone. This implies that intrinsic organization may reflect learned co-activation patterns. Relevant to interpreting why the manuscript finds intrinsic rather than extrinsic correlations between hands.


---

### 10. Representational models: a common framework for understanding encoding, pattern-component, and representational-similarity analysis

- **Authors**: Diedrichsen J, Kriegeskorte N

- **Year**: 2017

- **Journal**: PLoS Computational Biology, 13(4), e1005508

- **Link/DOI**: https://doi.org/10.1371/journal.pcbi.1005508

- **Relation score (1–10)**: 6

- **Verification**: Not PDF-verified


**Gap**

Multiple multivariate fMRI analysis methods (encoding models, RSA, pattern component models) existed but their relationships were unclear. A unifying framework was needed to understand when these methods agree or diverge.


**Task & Modality**

Theoretical/methodological paper with simulations and example fMRI data. Provides the analytical framework for decomposing activity patterns into components.


**Key results**

Encoding models, RSA, and pattern component models (PCM) are mathematically related and can be understood as different views of the same underlying second-moment matrix of activity patterns. This framework underpins the crossnobis dissimilarity and decomposition approach used in the manuscript. Shows how to properly test whether bimanual patterns can be decomposed into additive components.


---

### 11. Exploring interlimb constraints during bimanual graphic performance: effects of muscle grouping and direction

- **Authors**: Swinnen SP, Jardin K, Verschueren S, Meulenbroek R, Franz L, Dounskaia N, Walter CB

- **Year**: 1998

- **Journal**: Behavioural Brain Research, 90(1), 79-87

- **Link/DOI**: https://doi.org/10.1016/S0166-4328(97)00083-1

- **Relation score (1–10)**: 9

- **Verification**: Not PDF-verified


**Gap**

When people perform bimanual movements, interference can arise from similarity in muscle activation patterns (intrinsic) or similarity in movement direction (extrinsic). The relative contributions of these two constraint systems to bimanual interference were debated.


**Task & Modality**

Bimanual graphic (line-drawing) tasks where intrinsic (muscle) and extrinsic (spatial) congruency were independently manipulated. Behavioral kinematics.


**Key results**

Both intrinsic and extrinsic constraints affect bimanual coordination, but their relative dominance depends on task context. For rhythmic bimanual movements, intrinsic constraints (homologous muscle activation) were generally dominant, while extrinsic constraints became more important with visual feedback. This work established the intrinsic-extrinsic framework that the manuscript uses to test how contra- and ipsilateral activity patterns relate.


---

### 12. On the coordination of two-handed movements

- **Authors**: Kelso JAS, Southard DL, Goodman D

- **Year**: 1979

- **Journal**: Journal of Experimental Psychology: Human Perception and Performance, 5(2), 229-238

- **Link/DOI**: https://doi.org/10.1037/0096-1523.5.2.229

- **Relation score (1–10)**: 6

- **Verification**: Not PDF-verified


**Gap**

The fundamental principles governing bimanual coordination stability were not understood. Why do certain coordination patterns (in-phase/symmetric) come naturally while others (anti-phase) become unstable?


**Task & Modality**

Bimanual reaching movements with different amplitude and timing requirements. Behavioral analysis of temporal coupling between hands.


**Key results**

When the two hands had to move different distances, they tended to synchronize their movements temporally, with the shorter movement being slowed down. This seminal finding established that the motor system has a strong tendency toward temporal coupling between hands, consistent with the intrinsic coordination constraints that the manuscript's intrinsic neural correlations now provide a neural basis for.


---

### 13. Perceptual basis of bimanual coordination

- **Authors**: Mechsner F, Kerzel D, Knoblich G, Prinz W

- **Year**: 2001

- **Journal**: Nature, 414(6859), 69-73

- **Link/DOI**: https://doi.org/10.1038/35102060

- **Relation score (1–10)**: 7

- **Verification**: Not PDF-verified


**Gap**

Whether bimanual coordination is governed by perceptual/spatial (extrinsic) symmetry or by muscle-level (intrinsic) symmetry was debated. The dominant view favored motor/intrinsic constraints, but this study challenged that view.


**Task & Modality**

Bimanual circle-drawing with visual transformations that decoupled spatial from muscular symmetry. Behavioral coordination stability.


**Key results**

Argued that bimanual coordination is governed primarily by perceptual (extrinsic) symmetry, not muscle homology. When visual feedback was transformed so that non-homologous muscle patterns produced spatially symmetric output, coordination was stable. This directly contrasts with the manuscript's finding of intrinsic (not extrinsic) neural organization, suggesting the brain's internal coding may differ from what behavioral performance implies.


---

### 14. Muscle and movement representations in the primary motor cortex

- **Authors**: Kakei S, Hoffman DS, Strick PL

- **Year**: 1999

- **Journal**: Science, 285(5436), 2136-2139

- **Link/DOI**: https://doi.org/10.1126/science.285.5436.2136

- **Relation score (1–10)**: 9

- **Verification**: Not PDF-verified


**Gap**

Whether individual M1 neurons encode movement direction in extrinsic (spatial) or intrinsic (joint/muscle) coordinates was debated. Different wrist postures can dissociate these frames, revealing the true coordinate system.


**Task & Modality**

Monkey single-unit electrophysiology in M1 during wrist movements in different forearm postures (pronated vs. supinated) to dissociate extrinsic and intrinsic coordinates.


**Key results**

M1 neurons fall into distinct groups: some encode in extrinsic coordinates (direction tuning unchanged with posture), while others encode in intrinsic/muscle coordinates (tuning shifts with posture to follow muscle activation). This heterogeneity in M1 is directly relevant to the manuscript's question of whether bilateral representations are intrinsically or extrinsically organized. The wrist task is similar to the manuscript's.


---

### 15. Primary motor cortex is involved in bimanual coordination

- **Authors**: Donchin O, Gribova A, Steinberg O, Bergman H, Vaadia E

- **Year**: 1998

- **Journal**: Nature, 395(6699), 274-278

- **Link/DOI**: https://doi.org/10.1038/26220

- **Relation score (1–10)**: 8

- **Verification**: Not PDF-verified


**Gap**

During bimanual movements, it was unknown whether M1 neurons that respond to both hands do so in intrinsic (muscle) or extrinsic (spatial) coordinates. This is crucial for understanding bimanual interference at the neural level.


**Task & Modality**

Monkey single-unit electrophysiology in M1 during bimanual reaching. Compared neural tuning for ipsilateral movements in intrinsic vs. extrinsic coordinates.


**Key results**

Many M1 neurons showed significant modulation during ipsilateral movements. This landmark study established that M1 neurons carry bilateral information and are actively involved in bimanual coordination, not just contralateral hand control. Laid the groundwork for testing whether ipsilateral representations are organized in intrinsic or extrinsic coordinates at the population level, as the manuscript does with fMRI.


---

### 16. Integration of target and effector information in the human brain during reach planning

- **Authors**: Beurze SM, de Lange FP, Toni I, Medendorp WP

- **Year**: 2007

- **Journal**: Journal of Neurophysiology, 97(1), 188-199

- **Link/DOI**: https://doi.org/10.1152/jn.00456.2006

- **Relation score (1–10)**: 5

- **Verification**: Not PDF-verified


**Gap**

How parietal cortex integrates target location (extrinsic) with effector selection (left/right hand) during reaching was unclear. Do parietal areas use extrinsic, intrinsic, or intermediate reference frames?


**Task & Modality**

Reaching with left or right hand to visual targets during fMRI. Manipulated effector and target location to dissociate reference frames in posterior parietal cortex.


**Key results**

Posterior parietal cortex showed integration of target and effector information during reach planning. The representations showed a gradient from extrinsic to intrinsic coordinates across the parietal-to-premotor pathway. This is consistent with the manuscript's finding that parietal areas show bilateral representations and may explain why SPL areas show the intrinsic organization found.


---

### 17. Static and phasic cross-talk effects in discrete bimanual reversal movements

- **Authors**: Heuer H, Kleinsorge T, Spijkers W, Steglich C

- **Year**: 2001

- **Journal**: Journal of Motor Behavior, 33(1), 67-85

- **Link/DOI**: https://doi.org/10.1080/00222890109601904

- **Relation score (1–10)**: 7

- **Verification**: Not PDF-verified


**Gap**

Previous work treated intrinsic and extrinsic constraints as competing accounts, but their cooperative interaction was unexplored. Do these two reference frames combine additively or interact during bimanual movements?


**Task & Modality**

Bimanual reversal movements manipulating both intrinsic and extrinsic congruency simultaneously. Behavioral analysis of cross-talk and interference patterns.


**Key results**

Intrinsic and extrinsic constraints cooperate but with under-additivity, meaning the combined effect is less than the sum of individual effects. The intrinsic constraint was dominant for initial movement phases. This under-additive combination parallels the manuscript's sub-additive neural findings during bimanual movements, where ipsilateral representations diminish.


---

### 18. Ipsilateral motor cortex activity during unimanual hand movements relates to task complexity

- **Authors**: Verstynen T, Diedrichsen J, Albert N, Aparicio P, Ivry RB

- **Year**: 2005

- **Journal**: Journal of Neurophysiology, 93(3), 1209-1222

- **Link/DOI**: https://doi.org/10.1152/jn.00720.2004

- **Relation score (1–10)**: 8

- **Verification**: Not PDF-verified


**Gap**

The functional significance of ipsilateral M1 activity during unimanual movements was debated. Was it a genuine motor representation, an artifact of incomplete movement isolation, or a mirror of the contralateral pattern?


**Task & Modality**

Unimanual finger sequences of varying complexity during fMRI. Compared ipsilateral activation as a function of task complexity.


**Key results**

Confirmed genuine ipsilateral representations in M1 that scale with task complexity. Ipsilateral activation increased with more complex sequences. This validates the manuscript's finding of significant ipsilateral encoding in M1 and supports that it is not an artifact but reflects real motor information that scales with motor demands.


---

### 19. The cortical physiology of ipsilateral limb movements

- **Authors**: Bundy DT, Leuthardt EC

- **Year**: 2019

- **Journal**: Trends in Neurosciences, 42(11), 825-839

- **Link/DOI**: https://doi.org/10.1016/j.tins.2019.08.008

- **Relation score (1–10)**: 7

- **Verification**: Not PDF-verified


**Gap**

A comprehensive framework for understanding ipsilateral motor cortex activity was lacking. What functional roles does ipsilateral M1 serve, and how should we interpret ipsilateral activation in health and disease?


**Task & Modality**

Review article synthesizing fMRI, ECoG, TMS, and patient studies on ipsilateral motor cortex contributions to limb control.


**Key results**

Ipsilateral M1 activity reflects multiple sources: postural stabilization, bilateral coordination, and potential backup pathways. The review provides a framework for interpreting the manuscript's finding that ipsilateral encoding is present during unimanual movements but reduced during bimanual movements, suggesting context-dependent gating of ipsilateral representations.


---

### 20. The role of multiple contralesional motor areas for complex hand movements after internal capsular lesion

- **Authors**: Lotze M, Markert J, Sauseng P, Hoppe J, Plewnia C, Gerloff C

- **Year**: 2006

- **Journal**: Journal of Neuroscience, 26(22), 6096-6102

- **Link/DOI**: https://doi.org/10.1523/JNEUROSCI.4564-05.2006

- **Relation score (1–10)**: 4

- **Verification**: Not PDF-verified


**Gap**

After stroke damaging contralateral motor pathways, patients sometimes recover using ipsilateral pathways. Whether ipsilateral motor areas genuinely support movement recovery was clinically important and speaks to the latent bilateral representations found in healthy brains.


**Task & Modality**

Stroke patients with internal capsular lesions performed hand movements during fMRI. Compared ipsilateral motor area activation with healthy controls. Patient study.


**Key results**

Stroke patients showed dramatically increased contralesional (ipsilateral) motor cortex activation during movements of the affected hand. Multiple motor areas contributed, including M1, PMd, and SMA. Demonstrates that the ipsilateral representations documented in the manuscript (present but suppressed during bimanual movement) may serve as a latent backup system.


---

### 21. Effector-independent motor sequence representations exist in extrinsic and intrinsic reference frames

- **Authors**: Wiestler T, Waters-Metenier S, Diedrichsen J

- **Year**: 2014

- **Journal**: Journal of Neuroscience, 34(14), 5054-5064

- **Link/DOI**: https://doi.org/10.1523/JNEUROSCI.5363-13.2014

- **Relation score (1–10)**: 8

- **Verification**: Not PDF-verified


**Gap**

Whether ipsilateral movement representations in motor cortex are organized in intrinsic or extrinsic coordinates was unknown. If ipsilateral patterns share coordinate frames with contralateral patterns, they may reflect shared encoding principles.


**Task & Modality**

Finger sequence learning with transfer tests across hands, measured with fMRI. Multivariate analysis of ipsilateral and contralateral sequence representations.


**Key results**

Motor sequence representations exist in both extrinsic and intrinsic reference frames, and these persist across effectors. Learning a sequence with one hand created representations accessible by the other hand, organized in both coordinate systems. This dual coding is relevant to interpreting the manuscript's finding of intrinsic ipsilateral encoding.


---

### 22. Functional magnetic resonance imaging of motor cortex: hemispheric asymmetry and handedness

- **Authors**: Kim SG, Ashe J, Hendrich K, Ellermann JM, Merkle H, Ugurbil K, Georgopoulos AP

- **Year**: 1993

- **Journal**: Science, 261(5121), 615-617

- **Link/DOI**: https://doi.org/10.1126/science.8342027

- **Relation score (1–10)**: 6

- **Verification**: Not PDF-verified


**Gap**

Whether the degree of hemispheric lateralization in M1 varies with handedness and whether ipsilateral M1 activation is present during simple hand movements was unknown in the early fMRI era.


**Task & Modality**

Sequential finger-thumb opposition with dominant and non-dominant hands during fMRI. Compared contralateral vs. ipsilateral M1.


**Key results**

Found significant contralateral M1 activation for both hands, with additional ipsilateral M1 activation more prominent for the non-dominant hand. The ratio of contralateral to ipsilateral was about 3:1, remarkably similar to the 3.48x ratio reported in the manuscript. Established the baseline lateralization pattern that the manuscript builds upon.


---

### 23. Neural interactions between motor cortical hemispheres during bimanual and unimanual arm movements

- **Authors**: Cardoso de Oliveira S, Gribova A, Donchin O, Bergman H, Vaadia E

- **Year**: 2001

- **Journal**: European Journal of Neuroscience, 14(11), 1881-1896

- **Link/DOI**: https://doi.org/10.1046/j.0953-816x.2001.01801.x

- **Relation score (1–10)**: 8

- **Verification**: Not PDF-verified


**Gap**

The neural mechanism of bimanual interference in M1 was unknown. Do interhemispheric interactions between M1 populations change when the ipsilateral hand simultaneously moves, and could this cause interference?


**Task & Modality**

Monkey single-unit recording in bilateral M1 during unimanual and bimanual reaching. Analyzed interhemispheric neural correlations.


**Key results**

Interhemispheric correlations between M1 neurons changed significantly during bimanual compared to unimanual movements. This context-dependent modulation of inter-hemispheric communication could cause interference because the same population is influenced by the other hemisphere's activity. This parallels the manuscript's finding that ipsilateral encoding changes during bimanual movements.


---

### 24. Neuronal activity in cortical motor areas related to ipsilateral, contralateral, and bilateral digit movements of the monkey

- **Authors**: Tanji J, Okano K, Sato KC

- **Year**: 1988

- **Journal**: Journal of Neurophysiology, 60(1), 325-343

- **Link/DOI**: https://doi.org/10.1152/jn.1988.60.1.325

- **Relation score (1–10)**: 7

- **Verification**: Not PDF-verified


**Gap**

Whether single neurons in SMA and M1 are preferentially activated during bimanual vs. unimanual movements was unknown. SMA's specific role in bimanual coordination needed electrophysiological evidence.


**Task & Modality**

Monkey single-unit recording in SMA and M1 during ipsilateral, contralateral, and bilateral digit movements.


**Key results**

Found neurons in SMA specifically active during bimanual movements and not during either unimanual movement alone, providing direct evidence for bimanual-specific representations. M1 neurons were primarily contralateral-hand driven with less bimanual specificity. This SMA specialization is consistent with the manuscript's finding that the interaction component was present across regions including SMA.


---

### 25. Single-unit activity related to bimanual arm movements in the primary and supplementary motor cortices

- **Authors**: Donchin O, Gribova A, Steinberg O, Mitz AR, Bergman H, Vaadia E

- **Year**: 2002

- **Journal**: Journal of Neurophysiology, 88(6), 3498-3517

- **Link/DOI**: https://doi.org/10.1152/jn.00335.2001

- **Relation score (1–10)**: 8

- **Verification**: Not PDF-verified


**Gap**

Whether M1 neurons carry bimanual-specific information beyond what can be predicted from their unimanual responses was unknown. If bimanual activity is simply the sum of unimanual components, M1 may not actively coordinate bimanual movements.


**Task & Modality**

Monkey single-unit recording in M1 and SMA during unimanual and bimanual arm reaching. Compared actual bimanual responses to predictions from unimanual data.


**Key results**

Approximately 50% of M1 neurons showed bimanual responses that could not be predicted from their unimanual activities, demonstrating a genuine bimanual interaction at the single-neuron level. SMA showed even stronger bimanual-specific responses. This provides single-neuron evidence for the interaction component identified in the manuscript at the population level using fMRI.


---

### 26. Neural activity in primary motor and dorsal premotor cortex in reaching tasks with the contralateral versus ipsilateral arm

- **Authors**: Cisek P, Crammond DJ, Kalaska JF

- **Year**: 2003

- **Journal**: Journal of Neurophysiology, 89(2), 922-942

- **Link/DOI**: https://doi.org/10.1152/jn.00607.2002

- **Relation score (1–10)**: 7

- **Verification**: Not PDF-verified


**Gap**

The degree to which M1 and PMd neurons encode ipsilateral arm movements in primates was not systematically quantified. How lateralized are these areas, and does lateralization differ between M1 and premotor cortex?


**Task & Modality**

Monkey single-unit recording in M1 and PMd during unimanual reaching with contralateral vs. ipsilateral arm.


**Key results**

Both M1 and PMd showed significant modulation during ipsilateral arm movements, though weaker than contralateral. PMd showed more bilateral representation than M1, with similar directional tuning for both arms. The ipsilateral representations in PMd were consistent with intrinsic coordinates. This M1-PMd lateralization gradient matches the manuscript's fMRI finding.


---

### 27. Neural population dynamics during reaching

- **Authors**: Churchland MM, Cunningham JP, Kaufman MT, Foster JD, Nuyujukian P, Ryu SI, Shenoy KV

- **Year**: 2012

- **Journal**: Nature, 487(7405), 51-56

- **Link/DOI**: https://doi.org/10.1038/nature11129

- **Relation score (1–10)**: 3

- **Verification**: Not PDF-verified


**Gap**

Motor cortex activity during reaching was traditionally analyzed with static tuning curves, but a dynamics-based view was needed. How do M1 neural populations evolve over time during reaching?


**Task & Modality**

Monkey single-unit recording in M1 and PMd during delayed center-out reaching. Population dynamics analyzed with jPCA dimensionality reduction.


**Key results**

M1 population activity during reaching follows rotational dynamics consistent across conditions. This dynamical systems view provides context for understanding how the same M1 population might handle bilateral representations: the dynamics may occupy orthogonal subspaces for each hand, consistent with the manuscript's decomposition approach.


---

### 28. Motor cortex signals for each arm are mixed across hemispheres and neurons yet partitioned within the population response

- **Authors**: Ames KC, Churchland MM

- **Year**: 2019

- **Journal**: eLife, 8, e46159

- **Link/DOI**: https://doi.org/10.7554/eLife.46159

- **Relation score (1–10)**: 9

- **Verification**: Not PDF-verified


**Gap**

If M1 represents both hands, how does it avoid cross-talk? One hypothesis is that left and right hand information occupies orthogonal neural subspaces, preventing interference between the two hand representations.


**Task & Modality**

Monkey electrophysiology in M1 during unimanual and bimanual reaching. Dimensionality reduction to identify hand-specific subspaces.


**Key results**

Individual neurons were mixed in their selectivity for each arm, but at the population level, left and right arm representations occupied largely orthogonal neural subspaces. During bimanual movements, each hand's representation remained in its respective subspace, minimizing cross-talk. This orthogonality mechanism is directly relevant to the manuscript's discussion of how contralateral encoding remains stable while ipsilateral encoding disappears during bimanual movements.


---

### 29. Reorganization between preparatory and movement population responses in motor cortex

- **Authors**: Elsayed GF, Lara AH, Kaufman MT, Churchland MM, Cunningham JP

- **Year**: 2016

- **Journal**: Nature Communications, 7, 13239

- **Link/DOI**: https://doi.org/10.1038/ncomms13239

- **Relation score (1–10)**: 4

- **Verification**: Not PDF-verified


**Gap**

How motor cortex separates movement preparation from execution within the same neural population was unclear. The null-space hypothesis proposed that preparatory activity resides in a subspace that does not drive motor output.


**Task & Modality**

Monkey electrophysiology in M1 and PMd during delayed reaching. Subspace analysis separating preparatory from movement activity.


**Key results**

Preparatory and movement-related activity occupied orthogonal neural subspaces within the same population. This orthogonal subspace mechanism is the same principle that may allow M1 to represent both hands without interference. Relevant to understanding how the ipsilateral representation can exist without causing involuntary movements of the non-moving hand.


---

### 30. Supplementary motor area of the monkey's cerebral cortex: short- and long-term deficits after unilateral ablation and the effects of subsequent callosal section

- **Authors**: Brinkman C

- **Year**: 1984

- **Journal**: Journal of Neuroscience, 4(4), 918-929

- **Link/DOI**: https://doi.org/10.1523/JNEUROSCI.04-04-00918.1984

- **Relation score (1–10)**: 5

- **Verification**: Not PDF-verified


**Gap**

Whether SMA damage specifically impairs bimanual movements while sparing unimanual ones was unknown. Causal evidence for SMA's unique bimanual role was needed.


**Task & Modality**

Monkey lesion study: unilateral SMA ablation followed by testing unimanual and bimanual hand movements, then callosal section. Behavioral assessment.


**Key results**

SMA lesions caused severe deficits in bimanual coordination while unimanual movements remained relatively intact. Animals could still move each hand independently but could not coordinate them together. Subsequent callosal section worsened bimanual deficits. This provides causal evidence that SMA is necessary for bimanual coordination, supporting the manuscript's finding of bimanual-specific activity in SMA.


---

### 31. Clinical consequences of corticectomies involving the supplementary motor area in man

- **Authors**: Laplane D, Talairach J, Meininger V, Bancaud J, Orgogozo JM

- **Year**: 1977

- **Journal**: Journal of the Neurological Sciences, 34(3), 301-314

- **Link/DOI**: https://doi.org/10.1016/0022-510x(77)90148-4

- **Relation score (1–10)**: 5

- **Verification**: Not PDF-verified


**Gap**

Whether SMA damage in humans causes bimanual-specific deficits similar to monkey lesion studies was needed for cross-species validation of SMA's role.


**Task & Modality**

Neurological examination of patients with corticectomies involving the SMA. Assessment of unimanual vs. bimanual motor tasks. Patient study.


**Key results**

Patients with SMA lesions showed bimanual incoordination: normal unimanual movements but difficulty with bimanual tasks, especially asymmetric (non-mirror) movements. This human patient evidence complements monkey studies and is consistent with the manuscript's finding that SMA carries bimanual-specific interaction information.


---

### 32. The role of anterior cingulate cortex and precuneus in the coordination of motor behaviour

- **Authors**: Wenderoth N, Debaere F, Sunaert S, Swinnen SP

- **Year**: 2005

- **Journal**: European Journal of Neuroscience, 22(1), 235-246

- **Link/DOI**: https://doi.org/10.1111/j.1460-9568.2005.04176.x

- **Relation score (1–10)**: 7

- **Verification**: Not PDF-verified


**Gap**

Which cortical network supports bimanual coordination, and how does network activity change with coordination complexity? The specific roles of cingulate cortex and precuneus were unknown.


**Task & Modality**

Bimanual wrist movements with different coordination patterns (in-phase, anti-phase, multi-frequency ratios) during fMRI. Univariate activation mapping.


**Key results**

Bimanual coordination recruited anterior cingulate cortex, SMA, PMd, and precuneus beyond unimanual activation. More complex coordination patterns increased activation in these regions specifically. M1 activation was primarily contralateral. This network mapping is consistent with the manuscript's ROI selection and finding that premotor and medial wall areas show bimanual-specific activity.


---

### 33. Cerebellar and premotor function in bimanual coordination: parametric neural responses to spatiotemporal complexity and cycling frequency

- **Authors**: Debaere F, Wenderoth N, Sunaert S, van Hecke P, Swinnen SP

- **Year**: 2004

- **Journal**: NeuroImage, 21(4), 1416-1427

- **Link/DOI**: https://doi.org/10.1016/j.neuroimage.2003.12.011

- **Relation score (1–10)**: 7

- **Verification**: Not PDF-verified


**Gap**

The neural substrates distinguishing stable (in-phase) from unstable (anti-phase) bimanual coordination were not well characterized, particularly the cerebellar contribution.


**Task & Modality**

Bimanual cyclical wrist/forearm movements in in-phase, anti-phase, and 90-degree out-of-phase patterns during fMRI.


**Key results**

Anti-phase and complex coordination patterns produced greater activation in SMA, cerebellum, PMd, and secondary somatosensory cortex compared to in-phase. M1 showed similar activation for all patterns. This suggests the additional neural cost of non-intrinsic coordination is borne by premotor, supplementary, and cerebellar areas, consistent with the manuscript's interaction component findings.


---

### 34. Distribution of eye- and arm-movement-related neuronal activity in the SEF and in the SMA and pre-SMA of monkeys

- **Authors**: Fujii N, Mushiake H, Tanji J

- **Year**: 2002

- **Journal**: Journal of Neurophysiology, 87(4), 2158-2166

- **Link/DOI**: https://doi.org/10.1152/jn.00867.2001

- **Relation score (1–10)**: 3

- **Verification**: Not PDF-verified


**Gap**

How SMA and pre-SMA neurons are organized with respect to arm and eye movement control was unknown. Functional segregation within these areas could explain their roles in coordinating multi-effector and bimanual movements.


**Task & Modality**

Monkey single-unit recording in SMA, pre-SMA, and SEF during arm and eye movements. Analyzed spatial distribution of movement-related neurons.


**Key results**

SMA contained a high proportion of arm-movement-related neurons organized somatotopically, while pre-SMA contained more mixed arm/eye neurons. This functional organization supports SMA's role in limb coordination, including bimanual control, consistent with the manuscript's finding of bimanual-specific activity in medial motor areas.


---

### 35. Coordination of uncoupled bimanual movements by strictly timed interhemispheric connectivity

- **Authors**: Liuzzi G, Horniss V, Zimerman M, Gerloff C, Hummel FC

- **Year**: 2011

- **Journal**: Journal of Neuroscience, 31(25), 9111-9117

- **Link/DOI**: https://doi.org/10.1523/JNEUROSCI.0046-11.2011

- **Relation score (1–10)**: 8

- **Verification**: Not PDF-verified


**Gap**

The temporal dynamics of interhemispheric communication during bimanual coordination were unknown. Whether precisely timed connectivity between motor cortices is required for coordination needed causal testing.


**Task & Modality**

Paired-pulse TMS between motor cortices during bimanual finger movements. Measured temporal specificity of interhemispheric inhibition/facilitation.


**Key results**

Precisely timed interhemispheric connectivity was required for successful bimanual coordination. Disrupting this timing with TMS impaired coordination without affecting unimanual movements. This causal evidence supports the idea that the intrinsic coupling between hemispheric representations in the manuscript depends on callosal communication with precise timing.


---

### 36. Bimanual movement control: information processing and interaction effects

- **Authors**: Marteniuk RG, MacKenzie CL, Baba DM

- **Year**: 1984

- **Journal**: Quarterly Journal of Experimental Psychology, 36A(2), 335-365

- **Link/DOI**: https://doi.org/10.1080/14640748408402163

- **Relation score (1–10)**: 6

- **Verification**: Not PDF-verified


**Gap**

Whether bimanual interference (cross-talk) operates at the level of spatial plans or motor commands was unknown. Dissociating spatial from motoric interference reveals the locus of interlimb coupling.


**Task & Modality**

Bimanual aiming movements to targets of different sizes and distances. Behavioral analysis of movement time and accuracy for each hand.


**Key results**

When the two hands moved to targets of different difficulty, each hand's movement time was pulled toward the other. This assimilation effect demonstrated strong interlimb coupling at the motor planning level, consistent with shared representations. The behavioral interference matches the slight bimanual interference in the manuscript for incongruent conditions.


---

### 37. Callosotomy patients exhibit temporal uncoupling during continuous bimanual movements

- **Authors**: Kennerley SW, Diedrichsen J, Hazeltine E, Semjen A, Ivry RB

- **Year**: 2002

- **Journal**: Nature Neuroscience, 5(4), 376-381

- **Link/DOI**: https://doi.org/10.1038/nn822

- **Relation score (1–10)**: 8

- **Verification**: Not PDF-verified


**Gap**

Whether the corpus callosum is necessary for temporal coupling between hands during bimanual movements was debated. Some temporal coupling might be mediated subcortically.


**Task & Modality**

Callosotomy patients and controls performing continuous bimanual rhythmic movements. Measured temporal coupling. Patient study.


**Key results**

Callosotomy patients showed temporal uncoupling during continuous (but not discrete) bimanual movements, demonstrating that the corpus callosum is specifically necessary for ongoing temporal coordination. Discrete movements retained some subcortical coupling. The intrinsic coupling found in the manuscript at the cortical level likely depends on this callosal transmission.


---

### 38. Intermanual coordination: from behavioural principles to neural-network interactions

- **Authors**: Swinnen SP

- **Year**: 2002

- **Journal**: Nature Reviews Neuroscience, 3(5), 348-359

- **Link/DOI**: https://doi.org/10.1038/nrn807

- **Relation score (1–10)**: 7

- **Verification**: Not PDF-verified


**Gap**

A comprehensive review was needed to integrate behavioral principles of bimanual coordination with emerging neuroimaging and electrophysiology findings.


**Task & Modality**

Review covering behavioral experiments, fMRI, PET, lesion studies, and electrophysiology on bimanual coordination.


**Key results**

Identified key principles: bimanual coordination is constrained by both spatial and muscular symmetry; SMA and premotor areas are critical; the corpus callosum enables interhemispheric communication for asymmetric movements; M1 is primarily contralateral but carries ipsilateral information. Provides the conceptual framework within which the manuscript's findings should be interpreted.


---

### 39. Two hands, one brain: cognitive neuroscience of bimanual skill

- **Authors**: Swinnen SP, Wenderoth N

- **Year**: 2004

- **Journal**: Trends in Cognitive Sciences, 8(1), 18-25

- **Link/DOI**: https://doi.org/10.1016/j.tics.2003.10.017

- **Relation score (1–10)**: 6

- **Verification**: Not PDF-verified


**Gap**

How intermanual interactions during planning vs. execution contribute to bimanual interference needed synthesis across behavioral and imaging modalities.


**Task & Modality**

Review of behavioral and neuroimaging studies on bimanual skill learning and interference mechanisms.


**Key results**

Interference arises during both planning and execution, but the type differs. Planning-phase interference reflects shared representations in premotor/parietal areas (consistent with the manuscript's bilateral premotor encoding), while execution-phase interference involves M1 and subcortical structures. The intrinsic constraint is more dominant during execution, matching the intrinsic organization found.


---

### 40. Pattern component modeling: a flexible approach for understanding the representational structure of brain activity patterns

- **Authors**: Diedrichsen J, Yokoi A, Arbuckle SA

- **Year**: 2018

- **Journal**: NeuroImage, 180(Pt A), 119-133

- **Link/DOI**: https://doi.org/10.1016/j.neuroimage.2017.08.051

- **Relation score (1–10)**: 6

- **Verification**: Not PDF-verified


**Gap**

How to appropriately model and test the representational structure of brain activity patterns needed a flexible statistical framework beyond RSA.


**Task & Modality**

Methodological paper with finger movement fMRI examples. Introduces pattern component modeling (PCM) for testing representational models.


**Key results**

Introduced PCM as a likelihood-based framework for testing hypotheses about representational structure, complementing the crossnobis dissimilarity used in the manuscript. Cross-validated measures are essential to avoid bias, and component models can test specific hypotheses about how bimanual patterns decompose into additive vs. interactive components.


---

### 41. On the relations between the direction of two-dimensional arm movements and cell discharge in primate motor cortex

- **Authors**: Georgopoulos AP, Kalaska JF, Caminiti R, Massey JT

- **Year**: 1982

- **Journal**: Journal of Neuroscience, 2(11), 1527-1537

- **Link/DOI**: https://doi.org/10.1523/JNEUROSCI.02-11-01527.1982

- **Relation score (1–10)**: 3

- **Verification**: Not PDF-verified


**Gap**

How M1 neurons encode movement direction during reaching was a foundational question. Establishing directional tuning was needed for all subsequent work on contralateral and ipsilateral direction encoding.


**Task & Modality**

Monkey single-unit recording in M1 during center-out reaching to 8 targets. Characterized directional tuning curves and population vector decoding.


**Key results**

M1 neurons are broadly tuned for movement direction, each with a preferred direction. The population vector accurately predicted reach direction. This foundational work on directional encoding underpins the manuscript's analysis of directional encoding in motor regions. The same center-out paradigm is adapted in the manuscript for wrist movements.


---

### 42. Structure of population activity in primary motor cortex for single finger flexion and extension

- **Authors**: Arbuckle SA, Weiler J, Kirk EA, Rice CL, Schieber M, Pruszynski JA, Ejaz N, Diedrichsen J

- **Year**: 2020

- **Journal**: Journal of Neuroscience, 40(48), 9210-9223

- **Link/DOI**: https://doi.org/10.1523/JNEUROSCI.0999-20.2020

- **Relation score (1–10)**: 5

- **Verification**: PDF-verified


**Gap**

How the hand area of primary motor cortex (M1) is organized to control fine finger movements remains unclear, especially with respect to opposite movement directions. Because flexion and extension of the same finger never co-occur, a natural-statistics account would predict distinct representations—yet it was unknown whether population activity patterns in M1 separate flexion from extension, and how conclusions depend on measurement modality (fMRI vs spiking).


**Task & Modality**

Single-finger near-isometric flexion and extension task. Measured multivariate activity patterns with 7T fMRI in humans and compared them to neural spiking patterns from two monkeys performing the identical task; also related findings to muscular activity patterns.


**Key results**

Human fMRI patterns were distinct across different fingers, but were surprisingly similar for flexion versus extension of the same finger. In contrast, monkey M1 spiking patterns were distinct for both finger identity and movement direction (similar to muscular patterns). The fMRI–spiking discrepancy is consistent with an architecture where populations controlling flexion and extension of the same finger produce distinct outputs but share inputs and local recurrent processing and/or are spatially intermingled at a scale that makes opponent directions appear similar in coarse-scale fMRI. This provides testable hypotheses about how M1 organizes fine finger control.


---

### 43. Learned dynamics of reaching movements generalize from dominant to nondominant arm

- **Authors**: Criscimagna-Hemminger SE, Donchin O, Gazzaniga MS, Shadmehr R

- **Year**: 2003

- **Journal**: Journal of Neurophysiology, 89(1), 168-176

- **Link/DOI**: https://doi.org/10.1152/jn.00622.2002

- **Relation score (1–10)**: 4

- **Verification**: Not PDF-verified


**Gap**

Whether motor learning transfers between limbs and in which coordinate system was unknown. Transfer direction reveals the coordinate system of shared motor representations.


**Task & Modality**

Force-field adaptation with the dominant arm, then tested with the nondominant arm. Behavioral analysis of interlimb transfer.


**Key results**

Learning transferred from dominant to nondominant arm, and the pattern of transfer was consistent with generalization in extrinsic coordinates. This provided key evidence about how internal models are shared across limbs. Relevant to the manuscript's question of how representations relate between hands, whether in intrinsic or extrinsic coordinates.


---

### 44. Reach plans in eye-centered coordinates

- **Authors**: Batista AP, Buneo CA, Snyder LH, Andersen RA

- **Year**: 1999

- **Journal**: Science, 285(5425), 257-260

- **Link/DOI**: https://doi.org/10.1126/science.285.5425.257

- **Relation score (1–10)**: 2

- **Verification**: Not PDF-verified


**Gap**

In what coordinate frame does parietal cortex represent reach targets? Whether parietal reach plans are coded in eye-centered (extrinsic) or hand/body-centered (intrinsic) coordinates was debated.


**Task & Modality**

Monkey electrophysiology in parietal reach region during reaches with different eye and hand positions to dissociate reference frames.


**Key results**

Parietal reach region neurons encoded reach targets in eye-centered coordinates, an extrinsic frame of reference. This suggests that parietal cortex encodes goals in sensory coordinates before transformation into motor coordinates. Relevant to understanding the coordinate transformation gradient from extrinsic to intrinsic across the parietal-to-premotor pathway that the manuscript's intrinsic organization reflects.


---

### 45. Motor cortical representation of speed and direction during reaching

- **Authors**: Moran DW, Schwartz AB

- **Year**: 1999

- **Journal**: Journal of Neurophysiology, 82(5), 2676-2692

- **Link/DOI**: https://doi.org/10.1152/jn.1999.82.5.2676

- **Relation score (1–10)**: 2

- **Verification**: Not PDF-verified


**Gap**

Whether M1 neurons encode movement kinematics (direction, speed) or dynamics (force, torque) during reaching was debated. This distinction determines whether M1 operates in extrinsic or intrinsic coordinates.


**Task & Modality**

Monkey electrophysiology in M1 during continuous drawing movements. Compared neural correlates with kinematic (speed, direction) and dynamic (force) variables.


**Key results**

M1 neurons encoded both speed and direction of hand movement. The population showed a strong relationship to movement kinematics. This mixed representation means M1 does not purely operate in one coordinate frame, which helps explain why the manuscript finds intrinsic organization at the population level even though single neurons show mixed coding.


---
