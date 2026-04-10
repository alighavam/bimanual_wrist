# Bimanual coordination (applied fMRI/PET/DTI): annotated bibliography (~30)
Prepared: 2026-04-08

Format per paper:
- **Citation + link**
- **Gap (2 lines)**: what was unknown/contested and why this paper mattered.
- **Key points (3–4 lines)**: explanatory takeaways (not just what they did, but what it implies for interpreting bimanual fMRI).

---

## 1) Sadato et al., 1997, *J Neurosci* (PET)
**Sadato, N. et al.** “Role of the Supplementary Motor Area and the Right Premotor Cortex in the Coordination of Bimanual Finger Movements.” *Journal of Neuroscience* 17(24):9667–9674.  
Link: `https://doi.org/10.1523/JNEUROSCI.17-24-09667.1997`

- **Gap (2 lines)**: A core problem was separating “more movement” from “more coordination”: bimanual tasks typically co-vary with effort, attention, and bilateral M1 activity, making simple bimanual–unimanual subtraction misleading. This paper targets the *coordination-specific* component by contrasting mirror (default/easier) vs parallel (less natural/harder) bimanual patterns under matched pacing.
- **Key points (3–4 lines)**: The posterior **SMA** and **right PMd** show stronger activation for **parallel** compared to **mirror** bimanual movements, and also exceed unimanual activation in a second task—arguing these regions are not just “bilateral movement” markers. The interpretation is that SMA/PMd sit upstream of M1 as a coordination policy layer: they help implement non-default bimanual mappings and suppress mirror coupling. For modern analyses, this motivates focusing on medial-wall + premotor nodes (and their interhemispheric interactions) when your contrast is “independent hands” vs “coupled hands,” not simply “two hands” vs “one hand.”

---

## 2) Meyer-Lindenberg et al., 2002, *PNAS* (fMRI; coordination state transitions)
**Meyer-Lindenberg, A. et al.** “Transitions between Dynamical States of Differing Stability in the Human Brain.” *PNAS* 99:10948–10953.  
Link: `https://doi.org/10.1073/pnas.162114799`

- **Gap (2 lines)**: Coordination dynamics (anti-phase becoming unstable at higher frequency) were well described behaviorally, but the neural signature of a *state transition* was unclear: is it gradual load, or a distinct reconfiguration event? This paper frames the brain as moving between metastable coordination states and asks what changes when stability collapses.
- **Key points (3–4 lines)**: The important conceptual contribution is treating anti-phase→in-phase switching as a *qualitatively different brain event* rather than just “harder movement.” That matters because it legitimizes event-related analyses keyed to transition moments (not block averages), and it suggests that coordination failures may reflect network-level reweighting rather than local fatigue. In practice, it supports designing experiments that elicit transitions and then modeling them explicitly (e.g., time-locked network/ROI changes) rather than burying them as noise.

---

## 3) Koeneke et al., 2004, *NeuroImage* (uni vs bi; matched control tasks)
**Koeneke, S. et al.** “Bimanual versus unimanual coordination: what makes the difference?” *NeuroImage* 22(3):1336–1350.  
Link: `https://doi.org/10.1016/j.neuroimage.2004.03.012`

- **Gap (2 lines)**: Many papers implicitly assumed bimanual coordination recruits extra “bimanual-only” cortex, but comparisons were confounded by unmatched task structure. This paper asks, under carefully matched coordination demands, whether bimanual actions produce *additional* univariate activations beyond unimanual coordination.
- **Key points (3–4 lines)**: A central takeaway is negative but powerful: **bimanual>unimanual** did *not* robustly yield extra univariate activations when control tasks are properly constructed. That pushes interpretation away from “new regions appear” and toward “the same regions change how they interact or how efficiently they operate.” For your analyses, it implies that connectivity, multivariate pattern changes, or parametric scaling with coupling/independence may be more diagnostic than a blunt bimanual contrast.

---

## 4) Debaere et al., 2004, *Neuropsychologia* (learning a new bimanual pattern; fMRI)
**Debaere, F. et al.** “Changes in brain activation during the acquisition of a new bimanual coordination task.” *Neuropsychologia* 42(7):855–867.  
Link: `https://doi.org/10.1016/j.neuropsychologia.2003.12.010`

- **Gap (2 lines)**: It was unclear how the brain reallocates resources as a non-natural bimanual pattern becomes skilled: do “coordination regions” ramp up with learning, or do control networks step back as the task becomes automatized? This paper directly tracks learning-related activation changes rather than comparing novices vs experts across different contexts.
- **Key points (3–4 lines)**: Learning induces a characteristic *control-to-implementation* shift: regions associated with attentional/executive supervision and spatial control show decreasing engagement, while motor/subcortical/cerebellar circuitry shows increasing engagement as the skill stabilizes. Interpreting bimanual fMRI, this cautions that “more activation” is not always “better coordination”—early high activation may reflect compensation and control overhead. It also supports analyzing training as a trajectory and separating early acquisition from later automatization phases.

---

## 5) Liuzzi / Hummel lab, 2011, *J Neurosci* (effective connectivity; dpTMS + behavior)
**Liuzzi, G. et al.** “Coordination of Uncoupled Bimanual Movements by Strictly Timed Interhemispheric Connectivity.” *Journal of Neuroscience* 31(25):9111–9117.  
Link: `https://doi.org/10.1523/JNEUROSCI.0046-11.2011`

- **Gap (2 lines)**: Even if interhemispheric interactions matter, which pathways actually *predict* the ability to keep hands uncoupled (anti-phase, high frequency) was unresolved—especially with respect to timing (prep vs execution) and heterologous vs homologous connections. The paper ties physiological effective connectivity measures to specific coordination outcomes.
- **Key points (3–4 lines)**: The key explanatory point is that anti-phase success is linked to **fast, preparatory rPMd→lM1 facilitation**, not just generic “more coupling.” That suggests coordination failures can arise from mistimed interhemispheric gating, i.e., the signal arrives too late to stabilize the intended pattern. It also separates roles: homologous M1↔M1 interactions relate more to in-phase/unimanual regulation, while premotor-to-motor pathways contribute to independent control—useful when deciding which edges to test in connectivity models.

---

## 6) Johansen-Berg et al., 2007, *NeuroImage* (DTI; CC microstructure predicts bimanual skill)
**Johansen-Berg, H. et al.** “Integrity of white matter in the corpus callosum correlates with bimanual co-ordination skills.” *NeuroImage* 36(Suppl 2):T16–T21.  
Link: `https://doi.org/10.1016/j.neuroimage.2007.03.041`

- **Gap (2 lines)**: Behavioral work and lesion studies suggested the corpus callosum (CC) is important for bimanual coordination, but in healthy adults it was unclear whether *normal variation* in CC microstructure meaningfully predicts coordination ability. This paper tests whether structural connectivity differences explain why some people maintain asynchronous bimanual patterns better than others.
- **Key points (3–4 lines)**: The result is explanatory rather than merely correlational: the CC region whose FA predicts performance, when used as a tractography seed, connects to **SMA/cingulate motor** targets—i.e., plausible control hubs for bimanual coupling. That supports a model where the CC is not just “M1–M1 wiring,” but a route for medial-wall control signals that help decouple hands. For fMRI interpretation, it motivates treating CC integrity as a *moderator* of functional connectivity/behavior relationships.

---

## 7) Sisti et al., 2013, *Human Brain Mapping* (DTI; complexity + sensory feedback)
**Sisti, H.M. et al.** “Diffusion tensor imaging metrics of the corpus callosum in relation to bimanual coordination: Effect of task complexity and sensory feedback.” *Human Brain Mapping* 34(1):241–252.  
Link: `https://doi.org/10.1002/hbm.21429`

- **Gap (2 lines)**: Prior CC–bimanual links were typically shown for one task. It remained unclear whether CC microstructure matters most when coordination is complex (polyrhythms, non-1:1 ratios) and/or when feedback is limited (more internal control). This study manipulates both coordination complexity and augmented visual feedback to map when structure–behavior coupling is strongest.
- **Key points (3–4 lines)**: The main explanatory message is *conditional dependence*: CC microstructure relates to performance particularly when augmented feedback is removed, suggesting that interhemispheric pathways become more behaviorally limiting when you must rely on internal models/proprioception. It also supports functional parcellation of the CC: different subregions plausibly support different constraints (motor vs occipital/visual contributions). For designing bimanual fMRI, it implies sensory regime is not a nuisance variable; it changes which pathways constrain performance.

---

## 8) Sugawara et al., 2024/2025, *Experimental Brain Research* (7T fMRI; bimanual SRTT)
**Sugawara, S.K. et al.** “The left primary motor cortex and cerebellar vermis are critical hubs in bimanual sequential learning.” *Experimental Brain Research* 243(1):4.  
Link: `https://doi.org/10.1007/s00221-024-06944-2`

- **Gap (2 lines)**: Bimanual coordination learning blends at least two learning problems—learning stable bimanual “chords/postures” and learning sequence structure—but many studies collapse them into a single learning curve. This paper explicitly separates mirror vs parallel mappings and random vs sequence structure to ask what neural changes track chord formation vs sequence learning.
- **Key points (3–4 lines)**: The emphasis on **left M1** and **cerebellar vermis** as hubs is important because it argues coordination learning is not only premotor/SMA “planning”; primary motor and cerebellar circuits show systematic learning signatures. The connectivity finding (e.g., left M1 coupling with ACC) provides a mechanistic handle: learning changes not just regional activation but how control/monitoring systems integrate with motor output. For your work, it supports analyzing learning-related connectivity dynamics and testing whether parallel mapping selectively amplifies network reconfiguration.

---

## 9) Meister et al., 2010, *Cerebral Cortex* (handling temporally uncoupled bimanual movements)
**Meister, I.G. et al.** “How the Brain Handles Temporally Uncoupled Bimanual Movements.” *Cerebral Cortex* 20(12):2996–3004 (check exact pages in journal record).  
Link: `https://academic.oup.com/cercor/article/20/12/2996/369022`

- **Gap (2 lines)**: A recurring ambiguity in bimanual fMRI is whether difficulty reflects “more activation in the same motor regions” or a different configuration—especially when hands must be temporally uncoupled. This paper targets the uncoupled case to identify the broader control architecture needed to resist natural synchronization.
- **Key points (3–4 lines)**: The paper’s value is in emphasizing distributed contributions (premotor/parietal/prefrontal and interhemispheric interactions) for uncoupled timing, pushing beyond an M1-centric view. It supports interpreting anti-phase difficulty as a control problem—maintaining two partially independent controllers that still need coordinated goals. Methodologically, it encourages connectivity-aware analyses (e.g., PPI/SEM/DCM style reasoning) because the discriminative signal may lie in coupling patterns rather than peak activations.

---

## 10) Aramaki et al., 2006, *Cerebral Cortex* (spontaneous phase transition; event-related fMRI)
**Aramaki, Y. et al.** “Neural correlates of the spontaneous phase transition during bimanual coordination.” *Cerebral Cortex* 16(9):1338–1348.  
Link: `https://pubmed.ncbi.nlm.nih.gov/` (search by title; DOI not captured in current lookup)

- **Gap (2 lines)**: Many studies average across stable anti-phase blocks, but the most informative behavior may be the *moment* of breakdown and re-stabilization. This paper asks which brain areas are preferentially engaged at the transition itself versus those supporting steady-state execution of a chosen coordination mode.
- **Key points (3–4 lines)**: The transition-related recruitment of preSMA/rostral premotor/parietal regions suggests that switching is not a mere mechanical failure; it is accompanied by higher-level replanning or re-selection of a stable coordination solution. The dissociation from execution-related regions (SMA proper/caudal premotor/M1) provides a functional partition: one set maintains the motor pattern, another set handles instability and reconfiguration. Practically, it justifies modeling “transition epochs” as their own event type and testing whether right-hemisphere control asymmetry predicts who transitions sooner.

---

## 11) Rocca et al., 2007, *Human Brain Mapping* (in-phase vs anti-phase; posture/context manipulations)
**Rocca, M.A. et al.** “Influence of body segment position during in-phase and antiphase hand and foot movements: A kinematic and functional MRI study.” *Human Brain Mapping* 28:218–227.  
Link: `https://doi.org/10.1002/hbm.20271`

- **Gap (2 lines)**: Anti-phase is “harder,” but difficulty can also be induced by context (posture, spatial configuration), complicating simple interpretations of anti-phase>in-phase contrasts. This paper asks how coordination mode and body-segment position together shape the engaged networks, and whether the brain’s solution depends on spatial context.
- **Key points (3–4 lines)**: Anti-phase recruits broader sensorimotor and association networks, consistent with added control demands for maintaining an unstable pattern. The segment-position effect is important because it shows the control problem is partly representational: the mapping between intended coordination and body geometry changes the required computations. For bimanual fMRI design, it warns that posture/effector geometry can silently change “coordination load,” and those factors should be controlled or modeled.

---

## 12) Gerloff & Andres, 2002, *Acta Psychologica* (review; interhemispheric interaction framework)
**Gerloff, C. & Andres, F.G.** “Bimanual coordination and interhemispheric interaction.” *Acta Psychologica* 110(2–3).  
Link: `https://www.sciencedirect.com/science/article/abs/pii/S000169180200032X`

- **Gap (2 lines)**: Empirical results across imaging, EEG/MEG coherence, lesions, and stimulation were fragmented, making it hard to interpret any single fMRI finding as “the” coordination mechanism. This review articulates how interhemispheric coupling changes with task mode and learning stage, creating a functional vocabulary for interpreting bimanual tasks.
- **Key points (3–4 lines)**: The review’s explanatory value is that it positions bimanual coordination as an interhemispheric control problem with multiple possible coupling regimes (facilitation vs inhibition; homologous vs heterologous). It helps you interpret why some tasks show more SMA/PMd, others show more M1 or parietal, and why learning can invert activation patterns (control decreases, implementation increases). It also motivates multimodal designs (behavior + connectivity) because coordination may be “in the interactions,” not the activations.

---

## 13) Walsh et al., 2008, *NeuroImage* (network activation during bimanual movements)
**Walsh, R.R. et al.** “Network activation during bimanual movements in humans.” *NeuroImage* 43(3):540–553.  
Link: `https://pmc.ncbi.nlm.nih.gov/` (paper reported as PMCID: 2655207; DOI not captured in current lookup)

- **Gap (2 lines)**: Even when bimanual>unimanual univariate differences are subtle, it remains useful to characterize the *full network* reliably engaged by bimanual actions and how that network scales with task demands. This paper focuses on network-level activation patterns during bimanual movement to provide a reproducible map and baseline.
- **Key points (3–4 lines)**: The value here is giving a relatively complete catalog of regions consistently co-activated during bimanual movement (medial wall, premotor, parietal, cerebellar, subcortical). That helps anchor ROI selection and reduces the risk of over-interpreting one “hot” cluster as unique to your task. It also supports viewing bimanual behavior as coordination across distributed systems—timing, spatial mapping, monitoring—rather than a single locus.

---

## 14) Karabanov et al., 2020, *Frontiers in Human Neuroscience* (early bimanual skill learning; functional+structural plasticity)
**Karabanov, A. et al.** “Functional and Structural Plasticity Co-express in a Left Premotor Region During Early Bimanual Skill Learning.” *Frontiers in Human Neuroscience* 14:310.  
Link: `https://doi.org/10.3389/fnhum.2020.00310` (PMC: `https://pmc.ncbi.nlm.nih.gov/articles/PMC7456840/`)

- **Gap (2 lines)**: Many bimanual learning studies report functional activation changes without showing whether they reflect durable circuit remodeling or transient strategy shifts. This paper asks whether functional changes and structural (DTI) changes co-evolve in time during early bimanual skill acquisition, strengthening a plasticity-based interpretation.
- **Key points (3–4 lines)**: The co-expression of functional and structural change in a left premotor-related region supports the idea that early learning can rapidly re-tune communication pathways, not only reweight activations. It also provides a temporal claim: the biggest changes may occur early (scan 1→2), meaning “early learning” is a distinct regime worth isolating. For interpreting bimanual fMRI, it argues that premotor contributions can be both control-related and substrate-changing, depending on learning stage.

---

## 15) PLOS ONE, 2012 (motor imagery of bimanual everyday actions; connectivity emphasis)
**(PLOS ONE paper)** “Neural Activation and Functional Connectivity during Motor Imagery of Bimanual Everyday Actions.”  
Link: `https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0038506`

- **Gap (2 lines)**: Univariate activation contrasts often fail to cleanly isolate bimanual-specific effects, raising the question of whether bimanual coordination is primarily expressed as changes in *interregional coupling*. This paper uses motor imagery to reduce movement confounds and probes whether connectivity distinguishes bimanual coordination demands.
- **Key points (3–4 lines)**: The study illustrates a common pattern: limited or inconsistent bimanual-specific univariate activations, but clearer coordination-related changes in functional connectivity among parietal and premotor systems. That’s explanatory for why some fMRI bimanual studies “look null” unless connectivity is analyzed. It also suggests that even when overt movement is minimized (imagery), the coordination problem still recruits coupling mechanisms, supporting a cognitive-control framing.

---

## 16) Debaere et al., 2004, *NeuroImage* (parametric complexity & cycling frequency)
**Debaere, F. et al.** “Cerebellar and premotor function in bimanual coordination: parametric neural responses to spatiotemporal complexity and cycling frequency.” *NeuroImage* (2004).  
Link: `https://pubmed.ncbi.nlm.nih.gov/15050567/` (DOI reported in some sources as `10.1016/j.neuroimage.2003.12.011`; verify against journal page)

- **Gap (2 lines)**: “Complexity” and “speed” are often conflated: harder bimanual patterns are frequently performed at different frequencies, making it unclear whether SMA/PM/cerebellum effects reflect complexity per se or rate. This paper parametrically separates spatiotemporal complexity from cycling frequency to map which circuits scale with which demand dimension.
- **Key points (3–4 lines)**: The explanatory value is the dissociation: some regions scale with coordination complexity (reflecting planning/constraint handling), while others scale with frequency (reflecting execution/energetic demand). That helps interpret your own contrasts: if your manipulation changes frequency, you should expect certain “rate-sensitive” regions to respond even without increased coordination difficulty. Conversely, complexity-sensitive regions are stronger candidates for “coordination control” mechanisms and should align with behavioral coupling measures.

---

## 17) “Functional coupling of human cortical sensorimotor areas during bimanual skill acquisition” (Brain, 1999)
**(Brain, 1999 paper)** “Functional coupling of human cortical sensorimotor areas during bimanual skill acquisition.” *Brain* 122(5):855–870 (verify pages on journal record).  
Link: `https://academic.oup.com/brain/article/122/5/855/296636`

- **Gap (2 lines)**: Learning-related improvements are often described as regional activation increases/decreases, but coordination skill is plausibly about how areas cooperate. This paper targets *functional coupling* changes during skill acquisition, addressing whether learning is better described as network reconfiguration than local activation magnitude.
- **Key points (3–4 lines)**: The key explanatory point is that bimanual learning is accompanied by altered inter-area coupling within sensorimotor networks, consistent with the idea that skill is a coordination of controllers. This supports using coupling metrics (effective/functional connectivity, coherence analogs) as primary outcomes rather than secondary analyses. It also provides a conceptual bridge to later CC microstructure results: if callosal pathways constrain coupling, then learning-driven coupling changes should depend on those anatomical routes.

---

## 18) Aging/connectivity exemplar (Frontiers in Human Neuroscience, 2014)
**(Frontiers paper)** “Distant functional connectivity for bimanual finger coordination declines with aging: an fMRI and SEM exploration.”  
Link: `https://www.frontiersin.org/articles/10.3389/fnhum.2014.00251/full`

- **Gap (2 lines)**: Behavioral aging effects in bimanual tasks could come from weaker muscles, slower processing, or degraded coordination control; univariate fMRI alone may not disambiguate these. This paper uses connectivity modeling (SEM) to test whether aging primarily reduces long-range coordination links rather than simply reducing regional activation.
- **Key points (3–4 lines)**: The key explanatory contribution is that coordination deficits can be interpreted as *network disconnection*—reduced distant coupling—consistent with CC/white-matter vulnerability. Methodologically, it exemplifies how to connect behavior to directed connectivity hypotheses instead of relying on descriptive activation maps. It’s useful if you need an argument that bimanual coordination is especially sensitive to inter-area integration capacity.

---

## 19) Multiple sclerosis clinical anchor (J Neurosci, 2008)
**(J Neurosci, 2008 paper)** “Corpus Callosum and Bimanual Coordination in Multiple Sclerosis.”  
Link: `https://www.jneurosci.org/content/28/29/7248`

- **Gap (2 lines)**: Healthy-adult CC correlations are sometimes dismissed as “small effects.” A stronger mechanistic test is whether CC damage in disease systematically predicts bimanual impairment in ways consistent with interhemispheric coordination theories. This paper anchors the CC–bimanual link in a clinically meaningful disruption model.
- **Key points (3–4 lines)**: The explanatory value is convergent validity: when CC integrity is compromised, bimanual coordination suffers in predictable ways, supporting the CC as a bottleneck for interhemispheric timing/inhibition–facilitation control. It also helps interpret fMRI findings in patients: changes in activation may be compensatory responses to disconnection. For rehabilitation-oriented work, it motivates targeting interhemispheric communication and medial-wall motor areas rather than only ipsilesional M1 strengthening.

---

## 20) Bimanual coupling & left frontocentral network (2023; applied network association)
**(2023 paper; PubMed 37165733)** “Bimanual coupling is associated with left frontocentral network activity in a task-specific way.”  
Link: `https://pubmed.ncbi.nlm.nih.gov/37165733/`

- **Gap (2 lines)**: “Coupling” is a behavioral phenomenon observed across tasks, but it is often unclear whether it reflects a general trait-like neural signature or depends heavily on task demands. This paper asks whether coupling maps onto a consistent network—especially in left frontocentral circuitry—while allowing task dependence.
- **Key points (3–4 lines)**: The explanatory point is that coupling may not be a single mechanism; instead, left frontocentral involvement may rise specifically when task rules require suppressing default symmetry or maintaining a non-dominant mapping. That supports left-hemisphere dominance accounts in bimanual control and provides a plausible neural substrate for why some tasks show strong lateralization. For your analyses, it suggests relating coupling indices to network measures rather than expecting a universal SMA-only explanation.

---

## 21) Auditory-paced bimanual tapping lateralization exemplar (Clinical Neurophysiology abstract record)
**(Clinical Neurophysiology, 2016 abstract record)** “Lateralization of auditory-motor and somatosensory-motor loops during bimanual auditory-paced finger tapping.”  
Link: `https://www.sciencedirect.com/science/article/abs/pii/S1388245716301870`

- **Gap (2 lines)**: Many bimanual experiments pace movements, but pacing itself changes the control architecture (auditory–motor entrainment, sensory prediction). This record is useful for framing how pacing can induce lateralized loop engagement, potentially confounding interpretations of “bimanual coordination” effects.
- **Key points (3–4 lines)**: The key explanatory contribution is to treat pacing modality as a first-class factor: auditory-paced coordination may emphasize different circuits (auditory–motor coupling) than self-paced coordination (internal timing). That matters when comparing across studies or designing replication: different pacing could flip laterality or redistribute control burden. Use it to justify either standardizing pacing or explicitly modeling it when interpreting bimanual network activation.

---

## 22) Unilateral vs bilateral complex tapping activation exemplar (applied controls template)
**(PMC paper)** “Comparison of unilateral and bilateral complex finger tapping-related activation in premotor and primary motor cortex.”  
Link: `https://pmc.ncbi.nlm.nih.gov/articles/PMC6871138/`

- **Gap (2 lines)**: A persistent methodological gap is designing unimanual control tasks that truly match bilateral task complexity. This paper serves as a concrete example of how to compare unilateral and bilateral complex tapping while tracking premotor/M1 involvement under controlled conditions.
- **Key points (3–4 lines)**: The main value is methodological: it illustrates what changes are robust when bilaterality is introduced, and what changes disappear when complexity is matched—helping you calibrate expectations for bimanual>unimanual effects. It also provides practical ROI guidance (premotor vs M1) and a caution against over-attributing bilateral activation to “coordination” when it may reflect complexity or attention. If you need to defend your control condition design, this is a useful comparator.

---

## 23) Serrien et al., 2007, *Progress in Neurobiology* (action–cognition in coordination; review)
**Serrien, D.J., Ivry, R.B., & Swinnen, S.P.** (2007) *Progress in Neurobiology* review (see PDF).  
Link: `https://ivrylab.berkeley.edu/files/organized_pubs_pdfs/2007_serrien_ivry_swinnenprogr.pdf`

- **Gap (2 lines)**: Bimanual coordination studies often polarize into “motor-only” vs “cognitive control” accounts, but real tasks engage both; a coherent model was needed to explain why prefrontal/parietal engagement rises in complex coordination and learning. This review clarifies functional roles and reconciles evidence across modalities.
- **Key points (3–4 lines)**: The explanatory benefit is giving mechanistic language for interpreting why coordination difficulty recruits non-primary motor regions: response selection, conflict monitoring, and supervisory control become necessary when default coupling must be overridden. This helps interpret learning-related decreases in prefrontal activation as a reduction in control overhead rather than disengagement from the task. It also provides a theoretical rationale for connectivity analyses: cognitive-control regions may influence motor areas through timed gating signals, aligning with premotor→M1 facilitation findings.

---

## 24) “Distinct connectivity profiles predict different in-time processes of motor skill learning” (NeuroImage, 2021)
**(NeuroImage, 2021 paper)** “Distinct connectivity profiles predict different in-time processes of motor skill learning.”  
Link: `https://www.sciencedirect.com/science/article/pii/S1053811921005164`

- **Gap (2 lines)**: Skill learning is not unitary; early rapid improvements, mid-phase stabilization, and late refinement can have different neural signatures. This paper asks whether distinct connectivity profiles predict different temporal components of learning, offering a framework that can be applied to bimanual learning trajectories.
- **Key points (3–4 lines)**: The explanatory takeaway is that connectivity can act as a predictor of *which learning regime* a person is in, not only as a consequence of learning. Applied to bimanual tasks, this suggests modeling learning curves with multiple latent processes and mapping them onto different networks (premotor/parietal for rule implementation; M1/cerebellum for automatization). It helps you argue for richer learning models than pre/post contrasts.

---

## 25) “Brain–spinal cord interaction in long-term motor sequence learning” (NeuroImage, 2022)
**(NeuroImage, 2022 paper)** “Brain-spinal cord interaction in long-term motor sequence learning in human: An fMRI study.”  
Link: `https://www.sciencedirect.com/science/article/pii/S1053811922002397`

- **Gap (2 lines)**: Long-term learning might reorganize not only cortical networks but also descending control and spinal circuitry; many fMRI bimanual papers implicitly assume the locus is cortical. This paper expands the systems scope, arguing that sustained skill learning involves multi-level neural adaptation.
- **Key points (3–4 lines)**: The explanatory value is reminding you that stable bimanual skill can be supported by changes downstream of cortex, which can make cortical activation *decrease* even while performance improves. That matters when interpreting training effects: reduced cortical demand can reflect more efficient downstream implementation, not reduced engagement. It also suggests that “coordination difficulty” may partly reflect limitations in distributed control loops, not only cortical planning nodes.

---

## 26) Bimanual coordination constraints: Serrien, 2008, *Neuropsychologia* (behavioral/neural dynamics)
**Serrien, D.J.** “Coordination constraints during bimanual versus unimanual performance conditions.” *Neuropsychologia* 46:419–425.  
Link: `https://www.sciencedirect.com/science/article/abs/pii/S0028393207002953`

- **Gap (2 lines)**: Without a clear taxonomy of coordination constraints, imaging findings can be hard to interpret (is the subject failing due to spatial interference, temporal coupling, dominance effects, or feedback dependence?). This paper focuses on constraint differences between bimanual and unimanual contexts to sharpen the behavioral target of neuroimaging.
- **Key points (3–4 lines)**: The explanatory takeaway is that bimanual tasks introduce specific constraints (coupling tendencies, interference patterns, dominance asymmetries) that can be dissociated experimentally. This helps you map “what exactly is hard” in your bimanual paradigm to a hypothesized neural mechanism (e.g., conflict monitoring vs timing stabilization). It also motivates including explicit behavioral indices of coupling/phase error, because those indices determine which neural effect you should expect.

---

## 27) Frontiers (2019) pediatric lesion generalization: CC integrity vs bimanual coordination (CP)
**(Frontiers paper)** “Relationship Between Integrity of the Corpus Callosum and Bimanual Coordination in Children With Unilateral Spastic Cerebral Palsy.”  
Link: `https://www.frontiersin.org/articles/10.3389/fnhum.2019.00334/pdf`

- **Gap (2 lines)**: Adult healthy correlations do not guarantee causality or generality; developmental lesion contexts provide a stronger test of whether CC integrity is truly mechanistic for bimanual control. This paper asks whether CC integrity predicts bimanual coordination outcomes in a population where interhemispheric pathways are altered early in life.
- **Key points (3–4 lines)**: The explanatory contribution is that CC integrity remains behaviorally relevant under lesion constraints, supporting the idea that interhemispheric communication is not merely one of many redundant routes. It also provides clinically grounded bimanual measures that can be more ecologically valid than lab-only phase tasks. For interpreting fMRI in patient cohorts, it suggests anatomical constraints can dominate functional recruitment patterns, shaping both activation and connectivity.

---

## 28) TBI + CC microstructure + bimanual (young adults; diffusion study)
**(Monash record)** “Bimanual coordination and corpus callosum microstructure in young adults with traumatic brain injury: a diffusion tensor imaging study.”  
Link: `https://research.monash.edu/en/publications/bimanual-coordination-and-corpus-callosum-microstructure-in-young`

- **Gap (2 lines)**: Bimanual impairment after diffuse injury can arise from many sources; demonstrating a specific CC microstructure association supports an interhemispheric-communication mechanism. This study uses TBI as a natural experiment to test whether disrupted CC pathways map to coordination breakdown.
- **Key points (3–4 lines)**: The explanatory value is converging evidence across etiologies: different diseases/injuries that affect white matter also affect bimanual coordination, reinforcing that CC is a cross-condition bottleneck. It also helps interpret “why” fMRI patterns might differ in TBI: increased activation could reflect compensation for reduced coupling capacity. For applied work, it motivates combining diffusion + functional measures to separate structural limitation from strategy.

---

## 29) “Integrity of CC correlates with bimanual skill” extensions (task complexity + feedback; cross-paper synthesis)
**Synthesis entry (use alongside Johansen-Berg 2007 and Sisti 2013)**  
Links: `https://doi.org/10.1016/j.neuroimage.2007.03.041` and `https://doi.org/10.1002/hbm.21429`

- **Gap (2 lines)**: Individual papers often show one association, but the broader open question is *when* structure matters most: simple symmetric bimanual actions vs complex ratios, with vs without external feedback. This synthesis clarifies the boundary conditions so you can predict whether CC integrity should matter in your specific paradigm.
- **Key points (3–4 lines)**: Across these works, CC becomes most behaviorally predictive when coordination requires internal control and suppression of default coupling—precisely the conditions where interhemispheric timing and routing are stressed. That supports treating CC as a capacity constraint: if capacity is high, the system can maintain independence; if low, it collapses toward simpler coupled solutions. For fMRI, it implies that heterogeneous CC integrity across subjects can explain why some show strong premotor/SMA engagement or connectivity changes and others do not.

---

## 30) Behavioral anchor for interpreting phase-transition fMRI: Kelso (1984)
**Kelso, J.A.S.** (1984) classic coordination dynamics paper on anti-phase instability and transition to in-phase at higher frequencies.  
Link: (look up in your preferred database; not provided here)

- **Gap (2 lines)**: Imaging papers on bimanual phase transitions assume a behavioral law (anti-phase destabilizes with frequency) and need a canonical reference for that law. Without the behavioral anchor, neural “transition” findings can be misread as generic difficulty effects rather than a dynamical stability phenomenon.
- **Key points (3–4 lines)**: The explanatory point is that coordination patterns behave like attractor states with stability properties that change with control parameters (e.g., frequency). That provides a principled reason to design tasks that sweep frequency and to interpret sudden switches as state transitions rather than gradual fatigue. It also justifies modeling bimanual behavior with dynamical metrics (phase error, variability, switching probability) that can be related to neural markers.

