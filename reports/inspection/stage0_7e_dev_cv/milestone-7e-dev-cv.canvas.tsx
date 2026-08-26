import {
  BarChart,
  Callout,
  Card,
  CardBody,
  CardHeader,
  CollapsibleSection,
  Divider,
  Grid,
  H1,
  H2,
  H3,
  Pill,
  Row,
  Stack,
  Stat,
  Table,
  Text,
  useCanvasState,
} from "cursor/canvas";

type Tab = "results" | "schema" | "glossary" | "gaps";

const F1_CATEGORIES = [
  "Metadata only",
  "Multipath L1A",
  "M-value SGD-enet",
  "Gene-mean enet",
  "PCA-SVA",
  "M-value ridge",
  "Gene + direct",
  "Hier, no CpGPT",
  "M-value HGB",
  "Flat L1B",
];

const F1_VALUES = [0.659, 0.329, 0.324, 0.322, 0.29, 0.283, 0.27, 0.225, 0.103, 0.138];

const MAE_CATEGORIES = [
  "Metadata only",
  "M-value ridge",
  "M-value HGB",
  "Multipath L1A",
  "Gene + direct",
  "Gene-mean linear",
  "Hier, no CpGPT",
  "PCA-SVA",
  "Gene-mean enet",
];

const MAE_VALUES = [9.76, 10.77, 11.13, 11.49, 14.01, 15.59, 18.86, 19.03, 20.22];

export default function Milestone7eDevCv() {
  const [tab, setTab] = useCanvasState<Tab>("m7e-tab", "results");

  return (
    <Stack gap={24} style={{ padding: 24, maxWidth: 980 }}>
      <Stack gap={8}>
        <H1>Milestone 7E — architecture bake-off</H1>
        <Text tone="secondary">
          Frozen Age/Tissue/Sex Hub cohort: 13,548 samples, 327 studies, study-held-out
          3-fold CV. Winner rule: highest tissue macro-F1, then lowest age MAE in years.
        </Text>
      </Stack>

      <Grid columns={4} gap={16}>
        <Stat value="N-multipath-l1a" label="Selected architecture" tone="success" />
        <Stat value="0.329" label="Winner tissue macro-F1" />
        <Stat value="11.5 y" label="Winner age MAE" />
        <Stat value="0.659" label="Metadata-only F1 (ceiling)" tone="warning" />
      </Grid>

      <Callout tone="warning" title="Two-epoch neural budget">
        Deep Set models trained for 2 epochs on the first 8,192 CpG columns of a
        482,379-locus matrix. Linear and tree models on M-values used the same
        samples and the same locus prefix, but were not limited to 2 epochs. Do
        not read this as “trees beat Deep Sets.” It is an architecture-selection
        bake-off plus a check that methylation still carries signal under
        study-held-out splits.
      </Callout>

      <Row gap={8} wrap>
        <Pill active={tab === "results"} onClick={() => setTab("results")}>
          Results
        </Pill>
        <Pill active={tab === "schema"} onClick={() => setTab("schema")}>
          What went into each model
        </Pill>
        <Pill active={tab === "glossary"} onClick={() => setTab("glossary")}>
          Glossary
        </Pill>
        <Pill active={tab === "gaps"} onClick={() => setTab("gaps")}>
          Gaps and ROC
        </Pill>
      </Row>

      {tab === "results" && <ResultsPanel />}
      {tab === "schema" && <SchemaPanel />}
      {tab === "glossary" && <GlossaryPanel />}
      {tab === "gaps" && <GapsPanel />}

      <Divider />
      <Text tone="tertiary" size="small">
        Source: reports/inspection/stage0_7e_dev_cv/analysis.md · split
        hub-ats-7e-3fold-v1 · matrix-hub-age-tissue-sex-full-v1 · 26 Aug 2026
      </Text>
    </Stack>
  );
}

function ResultsPanel() {
  return (
    <Stack gap={20}>
      <H2>Tissue classification (higher is better)</H2>
      <Text tone="secondary" size="small">
        Macro-F1 averages F1 across tissue classes so rare tissues count as much as
        blood. Mean over 3 held-out study folds (neural arms also average 2 restarts).
      </Text>
      <BarChart
        categories={F1_CATEGORIES}
        series={[{ name: "Tissue macro-F1", data: F1_VALUES, tone: "info" }]}
        horizontal
        height={320}
      />

      <H2>Age prediction (lower is better)</H2>
      <Text tone="secondary" size="small">
        Mean absolute error in years. SGD elastic-net age is omitted: predictions
        diverged. Neural RMSE near 1 is standardized age, not years — MAE is the
        fair column.
      </Text>
      <BarChart
        categories={MAE_CATEGORIES}
        series={[{ name: "Age MAE (years)", data: MAE_VALUES, tone: "warning" }]}
        horizontal
        height={280}
      />

      <H2>Full comparison</H2>
      <Table
        headers={["Model", "Family", "Tissue F1", "Age MAE (y)", "Age R2"]}
        columnAlign={["left", "left", "right", "right", "right"]}
        striped
        rowTone={[
          "warning",
          "success",
          undefined,
          undefined,
          undefined,
          undefined,
          undefined,
          undefined,
          undefined,
          undefined,
        ]}
        rows={[
          ["Metadata (study + platform)", "control", "0.659", "9.76", "0.72"],
          ["Multipath late fusion (L1A)", "neural*", "0.329", "11.49", "0.62"],
          ["M-value SGD elastic-net", "classical", "0.324", "invalid", "—"],
          ["Gene-mean elastic-net", "transparent", "0.322", "20.22", "0.01"],
          ["M-value PCA-SVA + ridge", "classical", "0.290", "19.03", "-0.50"],
          ["M-value ridge + SGD-L2", "classical", "0.283", "10.77", "0.64"],
          ["Gene + direct CpG fusion", "neural*", "0.270", "14.01", "0.46"],
          ["Hierarchical gene, no CpGPT", "neural", "0.225", "18.86", "0.06"],
          ["M-value histogram boosting", "classical", "0.103", "11.13", "0.65"],
          ["Flat gene, Level-1", "neural", "0.138", "18.38", "-0.20"],
        ]}
      />
      <Text tone="tertiary" size="small">
        neural* = reported numbers are late-fusion linear heads on region means,
        not saved Deep Set MBS matrices. Metadata is a leakage alarm, not a model
        to ship.
      </Text>
    </Stack>
  );
}

function SchemaPanel() {
  return (
    <Stack gap={16}>
      <H2>Data and targets</H2>
      <Table
        headers={["Item", "Definition"]}
        rows={[
          ["ATS cohort", "Age / tissue / sex Hub GSM-union, 13,548 samples"],
          ["Matrix", "matrix-hub-age-tissue-sex-full-v1, GRCh38, Illumina arrays"],
          ["Split", "hub-ats-7e-3fold-v1 — same study never in train and test"],
          ["Age", "Continuous years. Metric: MAE in years"],
          ["Tissue", "Many classes. Metric: macro-F1, balanced accuracy"],
          ["Sex", "Two classes. Metric: AUROC. ROC applies here"],
          ["Unknown labels", "Missing is unknown, never treated as a control"],
        ]}
      />

      <H2>What each arm combines</H2>
      <Table
        headers={["Arm", "Inputs", "Pooling", "Head"]}
        rows={[
          ["T-mean / T-enet", "Mean beta per gene", "Mean", "Ridge or elastic-net + logistic"],
          ["N-flat-gene", "beta, M, optional z, optional CpGPT", "Max Deep Set", "Joint linear heads"],
          ["N-hier-gene", "Same + region type + residual CpGs", "Two-level Deep Set", "Joint linear heads"],
          ["N-rbs / N-tbs", "Same, restricted to CGI or tiles", "Flat Deep Set", "Joint linear heads"],
          ["N-gene-direct", "Gene means + per-CpG elastic-net preds", "Late linear", "Ridge / logistic"],
          ["N-multipath", "Gene + RBS + TBS means + direct", "Late linear", "Ridge / logistic"],
          ["C-mvalue-*", "M = log2(beta/(1-beta)), 8,192 loci", "None or trees", "Ridge / SGD / HGB"],
          ["C-mvalue-sva", "M residualized on 10 train PCs", "None", "Ridge / SGD logistic"],
          ["C-metadata", "Study id + platform only", "—", "Ridge / logistic"],
        ]}
      />

      <H3>Encoder (flat and hierarchical, matched)</H3>
      <Text>
        GELU, dropout 0.1, LayerNorm, CpG hidden size 64. Level-1 A = raw M; B =
        fold-fitted robust z of M. CpGPT is a frozen locus embedding, not
        sample methylation.
      </Text>
    </Stack>
  );
}

function GlossaryPanel() {
  return (
    <Stack gap={8}>
      <H2>Terms</H2>
      <CollapsibleSection title="MBS / RBS / TBS" defaultOpen>
        <Text>
          Methylation, Regulatory, and Tile Burden Scores. One number (or vector)
          per gene, CGI/regulatory region, or intergenic 50-CpG tile. Tiles are
          not assigned to the nearest gene.
        </Text>
      </CollapsibleSection>
      <CollapsibleSection title="Beta vs M-value">
        <Text>
          Beta is the methylated fraction (about 0-1). M-value is
          log2(beta / (1 - beta)), closer to Gaussian, standard for linear models.
        </Text>
      </CollapsibleSection>
      <CollapsibleSection title="Deep set / flat / hierarchical">
        <Text>
          A Deep Set pools a variable-length list of CpGs (order does not matter).
          Flat: CpGs to gene. Hierarchical: CpGs to region to gene, plus leftover
          CpGs.
        </Text>
      </CollapsibleSection>
      <CollapsibleSection title="Level-1, MAD, CpGPT">
        <Text>
          Level-1 is a train-fold-only robust z-score of M using median and MAD
          (median absolute deviation). CpGPT is a frozen DNA-language vector for
          each locus, shared across samples.
        </Text>
      </CollapsibleSection>
      <CollapsibleSection title="Late fusion vs direct">
        <Text>
          Late fusion trains branches separately then glues features with a linear
          head. Direct is a per-CpG penalised model whose predictions are extra
          columns — no gene pooling.
        </Text>
      </CollapsibleSection>
      <CollapsibleSection title="Metrics: macro-F1, MAE, AUROC">
        <Text>
          Macro-F1: mean F1 across classes. MAE: mean |predicted − true| years.
          AUROC: ranking quality; 0.5 is chance, 1 is perfect. ROC is for sex and
          one-vs-rest tissue, not for age.
        </Text>
      </CollapsibleSection>
      <CollapsibleSection title="SVA, HGB, elastic-net, ridge">
        <Text>
          SVA here: PCA on train M-values, residualize, then linear models. HGB:
          sklearn histogram gradient boosting (LightGBM family; LightGBM was not
          added as a dependency). Elastic-net: L1+L2. Ridge: L2 only. Linear
          classifiers used SGD because SAGA/coordinate descent did not finish on
          8,192 loci.
        </Text>
      </CollapsibleSection>
    </Stack>
  );
}

function GapsPanel() {
  return (
    <Stack gap={16}>
      <H2>ROC (M-value HGB, fold 0 held-out studies)</H2>
      <Text tone="secondary">
        One-versus-rest: one tissue vs all others. Age has no ROC. Three of five
        plotted tissues sit on the chance line in this fold (absent or unlearned
        under the 40-iteration tree budget).
      </Text>
      <Table
        headers={["Target", "AUROC", "How to read"]}
        columnAlign={["left", "right", "left"]}
        rows={[
          ["Sex", "0.75", "Binary ranking of female vs male"],
          ["Whole blood vs rest", "0.84", "HGB can rank a common tissue"],
          ["Breast vs rest", "0.76", "Better than chance"],
          ["Kidney / lung / semen vs rest", "0.50", "Chance in this holdout"],
        ]}
      />

      <H2>What is still missing</H2>
      <Card>
        <CardHeader>Evaluation gaps, not a crashed trainer</CardHeader>
        <CardBody>
          <Stack gap={8}>
            <Text>
              The 90-cell neural bake-off finished (3 folds x arms). These are
              quality gaps:
            </Text>
            <Text>
              1. Neural nets: 2 epochs, 8,192-locus prefix — not a converged Deep
              Set.
            </Text>
            <Text>
              2. Multipath numbers are linear fusion of region means, not neural
              MBS matrices.
            </Text>
            <Text>
              3. Full 482,379 CpGs were not used (memory). Comparison is matched
              to the neural prefix.
            </Text>
            <Text>
              4. No Bioconductor sva; 10 train PCs. No LightGBM package; HGB is
              the same tree family and was under-trained for 47 tissues (F1 0.10)
              while age MAE 11.1 years is competitive.
            </Text>
            <Text>
              5. SGD elastic-net age diverged. Tissue logistic from that family
              is usable (F1 0.324).
            </Text>
          </Stack>
        </CardBody>
      </Card>

      <Callout tone="success" title="Recommendation for Milestone 7">
        Keep multi-path topology (gene + RBS + TBS + direct) with late fusion.
        Revisit joint DeepRVAT heads after a longer train on these same frozen
        folds. Do not ship metadata-only. Do not start 5x6 OOF until unknown
        disease/cancer labels stay unknown (7E' is done).
      </Callout>
    </Stack>
  );
}
