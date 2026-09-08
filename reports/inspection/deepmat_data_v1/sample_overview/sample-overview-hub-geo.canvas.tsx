import {
  BarChart,
  Callout,
  Card,
  CardBody,
  CardHeader,
  Grid,
  H1,
  H2,
  Pill,
  Row,
  Stack,
  Stat,
  Table,
  Text,
  useCanvasState,
} from "cursor/canvas";

type Tab = "overview" | "phenotypes" | "tissues" | "diseases";

/** Aggregates from reports/inspection/deepmat_data_v1/sample_overview/summary.json */
const N_SAMPLES = 173076;
const LANE = { hub: 34234, ewas: 169176, both: 30334, geo: 90751 };
const PACKS = [
  { pack: "disease", n: 12218 },
  { pack: "cancer", n: 10101 },
  { pack: "age", n: 8374 },
  { pack: "tissue", n: 5323 },
  { pack: "blood", n: 3402 },
  { pack: "sex", n: 2978 },
  { pack: "bmi", n: 2070 },
  { pack: "brain", n: 1997 },
  { pack: "ancestry", n: 1380 },
];
const COVERAGE = [
  { field: "tissue", n: 167756 },
  { field: "pubmed_ids", n: 159483 },
  { field: "sex", n: 82833 },
  { field: "age_years", n: 58054 },
  { field: "cell_component", n: 3280 },
  { field: "bmi", n: 2070 },
  { field: "brain_label", n: 1997 },
  { field: "ancestry_label", n: 1380 },
];
const PLATFORMS = [
  { platform: "HM450", n: 109903 },
  { platform: "EPIC", n: 56711 },
  { platform: "EPICv2", n: 6415 },
  { platform: "unknown", n: 47 },
];
const SPECIES = [
  { status: "human", n: 172721 },
  { status: "non_human", n: 355 },
];
const SPECIES_SOURCE = [
  { source: "geo_soft", n: 90751 },
  { source: "geo_series_taxon", n: 68567 },
  { source: "project_namespace_human", n: 6314 },
  { source: "hub_baseline_assumed_human", n: 6269 },
  { source: "ewas_datahub_assumed_human", n: 1175 },
];
const PLATFORM_SOURCE = [
  { source: "explicit", n: 167357 },
  { source: "ewas_db_n_probes", n: 3527 },
  { source: "datahub_metadata_json", n: 2142 },
  { source: "sample_id_935k_hint", n: 3 },
  { source: "unresolved", n: 47 },
];
const SEX = [
  { label: "Female", n: 36151 },
  { label: "Male", n: 34189 },
  { label: "M", n: 6581 },
  { label: "F", n: 5904 },
];
const ANCESTRY = [
  { label: "African American", n: 237 },
  { label: "Hispanic - Mexican", n: 227 },
  { label: "Hispanic", n: 168 },
  { label: "Chinese", n: 167 },
  { label: "white", n: 159 },
  { label: "Caucasian - European", n: 144 },
  { label: "East Asian", n: 74 },
  { label: "Northern European", n: 73 },
  { label: "African", n: 72 },
];
const TISSUES = [
  { tissue: "whole blood", n: 35741 },
  { tissue: "brain - tumor", n: 11452 },
  { tissue: "leukocyte", n: 7456 },
  { tissue: "peripheral blood mononuclear cell", n: 7440 },
  { tissue: "brain", n: 3828 },
  { tissue: "breast", n: 3519 },
  { tissue: "lung", n: 3108 },
  { tissue: "cord blood", n: 2882 },
  { tissue: "liver", n: 2470 },
  { tissue: "semen", n: 2100 },
  { tissue: "kidney", n: 2056 },
  { tissue: "saliva", n: 1798 },
  { tissue: "placenta", n: 1721 },
  { tissue: "CD14+ monocyte", n: 1624 },
  { tissue: "bone marrow", n: 1544 },
];
const BRAIN = [
  { region: "cerebellum", n: 300 },
  { region: "dlPFC", n: 245 },
  { region: "superior temporal gyrus", n: 179 },
  { region: "frontal lobe", n: 155 },
  { region: "frontal cortex", n: 133 },
  { region: "occipital lobe", n: 124 },
  { region: "temporal lobe", n: 121 },
  { region: "entorhinal cortex", n: 101 },
];

/** Hub disease pack (nine-pack); diagnoses from case phenotype_value */
const DISEASE = {
  nNinePack: 34234,
  nPack: 12218,
  nCase: 5264,
  nControl: 6930,
  nAdjacent: 24,
  nUsable: 12194,
  nLabels: 28,
  fracOfPack: 0.4308,
  fracOfNinePack: 0.1538,
};
const CANCER = {
  nPack: 10101,
  nCase: 7157,
  nControl: 1920,
  nAdjacent: 1024,
  fracOfPack: 0.7085,
};
const DISEASES = [
  { disease: "Alzheimer's disease", n: 945, pct: 17.95 },
  { disease: "schizophrenia", n: 536, pct: 10.18 },
  { disease: "systemic lupus erythematosus", n: 341, pct: 6.48 },
  { disease: "Parkinson's disease", n: 333, pct: 6.33 },
  { disease: "ulcerative colitis", n: 258, pct: 4.9 },
  { disease: "multiple sclerosis", n: 228, pct: 4.33 },
  { disease: "rheumatoid arthritis", n: 225, pct: 4.27 },
  { disease: "psoriasis", n: 211, pct: 4.01 },
  { disease: "stroke", n: 204, pct: 3.88 },
  { disease: "childhood asthma", n: 202, pct: 3.84 },
  { disease: "intellectual disability and congenital anomalies", n: 200, pct: 3.8 },
  { disease: "Crohn's disease", n: 197, pct: 3.74 },
  { disease: "asthma", n: 194, pct: 3.69 },
  { disease: "preeclampsia", n: 179, pct: 3.4 },
  { disease: "Huntington's disease", n: 170, pct: 3.23 },
  { disease: "systemic insulin resistance", n: 115, pct: 2.18 },
  { disease: "respiratory allergy", n: 104, pct: 1.98 },
  { disease: "type 2 diabetes", n: 92, pct: 1.75 },
  { disease: "panic disorder", n: 89, pct: 1.69 },
  { disease: "Silver Russell syndrome", n: 79, pct: 1.5 },
  { disease: "Graves' disease", n: 73, pct: 1.39 },
  { disease: "spina bifida", n: 59, pct: 1.12 },
  { disease: "Down syndrome", n: 54, pct: 1.03 },
  { disease: "autism spectrum disorder", n: 48, pct: 0.91 },
  { disease: "Sjogren's syndrome", n: 48, pct: 0.91 },
  { disease: "Kabuki syndrome", n: 40, pct: 0.76 },
  { disease: "nephrogenic rest", n: 22, pct: 0.42 },
  { disease: "systemic sclerosis", n: 18, pct: 0.34 },
];

export default function SampleOverviewHubGeo() {
  const [tab, setTab] = useCanvasState<Tab>("overview-tab", "overview");

  return (
    <Stack gap={24} style={{ padding: 24, maxWidth: 980 }}>
      <Stack gap={8}>
        <H1>Sample overview — Hub + EWAS_db</H1>
        <Text tone="secondary">
          Full catalog census (n={N_SAMPLES.toLocaleString()}). GSM = sample, GSE
          = study, GPL = platform. Project copy + PDF:
          reports/inspection/deepmat_data_v1/sample_overview/
        </Text>
      </Stack>

      <Row gap={8} wrap>
        <Pill active={tab === "overview"} onClick={() => setTab("overview")}>
          Overview
        </Pill>
        <Pill
          active={tab === "phenotypes"}
          onClick={() => setTab("phenotypes")}
        >
          Phenotypes
        </Pill>
        <Pill active={tab === "tissues"} onClick={() => setTab("tissues")}>
          Tissues
        </Pill>
        <Pill active={tab === "diseases"} onClick={() => setTab("diseases")}>
          Diseases
        </Pill>
      </Row>

      {tab === "overview" && <OverviewPanel />}
      {tab === "phenotypes" && <PhenotypesPanel />}
      {tab === "tissues" && <TissuesPanel />}
      {tab === "diseases" && <DiseasesPanel />}
    </Stack>
  );
}

function OverviewPanel() {
  return (
    <Stack gap={20}>
      <Callout tone="info">
        Species filled for every row (GEO SOFT / series taxon / Hub & project
        priors). Unknown platforms inferred from DataHub metadata, GPL maps, and
        EWAS_db probe counts (
        {PLATFORM_SOURCE.find((p) => p.source === "ewas_db_n_probes")?.n}{" "}
        samples). {PLATFORMS.find((p) => p.platform === "unknown")?.n}{" "}
        custom/truncated files remain unresolved.
      </Callout>

      <Grid columns={4} gap={12}>
        <Stat value={N_SAMPLES.toLocaleString()} label="Catalog samples" />
        <Stat value={LANE.hub.toLocaleString()} label="Hub baseline" />
        <Stat value={SPECIES[0].n.toLocaleString()} label="Human" />
        <Stat value={SPECIES[1].n.toLocaleString()} label="Non-human" />
      </Grid>

      <Card>
        <CardHeader>Lane membership</CardHeader>
        <CardBody>
          <BarChart
            categories={["Hub baseline", "EWAS_db", "Hub ∩ EWAS", "GEO SOFT"]}
            series={[
              {
                name: "Samples",
                data: [LANE.hub, LANE.ewas, LANE.both, LANE.geo],
              },
            ]}
            height={220}
          />
        </CardBody>
      </Card>

      <Grid columns={2} gap={16}>
        <Card>
          <CardHeader>Platform (after inference)</CardHeader>
          <CardBody>
            <BarChart
              categories={PLATFORMS.map((p) => p.platform)}
              series={[{ name: "Samples", data: PLATFORMS.map((p) => p.n) }]}
              height={220}
            />
            <Table
              headers={["Source", "N"]}
              rows={PLATFORM_SOURCE.map((p) => [
                p.source,
                p.n.toLocaleString(),
              ])}
            />
          </CardBody>
        </Card>
        <Card>
          <CardHeader>Species provenance</CardHeader>
          <CardBody>
            <BarChart
              categories={SPECIES.map((s) => s.status)}
              series={[{ name: "Samples", data: SPECIES.map((s) => s.n) }]}
              height={160}
            />
            <Table
              headers={["Source", "N"]}
              rows={SPECIES_SOURCE.map((s) => [
                s.source,
                s.n.toLocaleString(),
              ])}
            />
          </CardBody>
        </Card>
      </Grid>

      <Stack gap={8}>
        <H2>Disease overview (Hub pack)</H2>
        <Text tone="secondary">
          {DISEASE.nLabels} diagnoses ·{" "}
          {(100 * DISEASE.fracOfPack).toFixed(1)}% of pack are patients ·{" "}
          {(100 * DISEASE.fracOfNinePack).toFixed(1)}% of nine-pack. Full list
          on the Diseases tab.
        </Text>
      </Stack>
      <Grid columns={4} gap={12}>
        <Stat value={DISEASE.nPack.toLocaleString()} label="Disease pack" />
        <Stat value={DISEASE.nCase.toLocaleString()} label="Cases" />
        <Stat value={DISEASE.nControl.toLocaleString()} label="Controls" />
        <Stat
          value={`${(100 * DISEASE.fracOfPack).toFixed(1)}%`}
          label="Patient fraction"
        />
      </Grid>
      <Card>
        <CardHeader>Top diagnoses (cases)</CardHeader>
        <CardBody>
          <BarChart
            horizontal
            categories={DISEASES.slice(0, 10).map((d) => d.disease)}
            series={[
              {
                name: "Case samples",
                data: DISEASES.slice(0, 10).map((d) => d.n),
              },
            ]}
            height={320}
          />
        </CardBody>
      </Card>

      <Card>
        <CardHeader>Hub pack membership</CardHeader>
        <CardBody>
          <BarChart
            categories={PACKS.map((p) => p.pack)}
            series={[{ name: "Samples", data: PACKS.map((p) => p.n) }]}
            height={240}
          />
        </CardBody>
      </Card>
    </Stack>
  );
}

function PhenotypesPanel() {
  return (
    <Stack gap={20}>
      <H2>Phenotype coverage</H2>
      <Text tone="secondary">
        Non-null counts after Hub-wins / GEO-fill. BMI and ancestry are real Hub
        labels (freeze-and-reuse heads wired); brain is region catalogue only.
      </Text>
      <Card>
        <CardHeader>Coverage bars</CardHeader>
        <CardBody>
          <BarChart
            categories={COVERAGE.map((c) => c.field)}
            series={[{ name: "Samples", data: COVERAGE.map((c) => c.n) }]}
            height={280}
          />
        </CardBody>
      </Card>
      <Grid columns={2} gap={16}>
        <Stack gap={8}>
          <H2>Sex labels</H2>
          <Table
            headers={["Label", "N"]}
            rows={SEX.map((s) => [s.label, s.n.toLocaleString()])}
          />
        </Stack>
        <Stack gap={8}>
          <H2>Ancestry (Hub)</H2>
          <Table
            headers={["Label", "N"]}
            rows={ANCESTRY.map((a) => [a.label, a.n.toLocaleString()])}
          />
        </Stack>
      </Grid>
      <Callout tone="neutral">
        Sex still has raw F/M vs Female/Male strings — ontology collapse is a
        separate hygiene task. Age coverage is lower than tissue because many
        GEO-only rows lack numeric age.
      </Callout>
    </Stack>
  );
}

function TissuesPanel() {
  return (
    <Stack gap={20}>
      <H2>Tissues overview</H2>
      <Text tone="secondary">
        Harmonized tissue labels across Hub + GEO (Hub wins). Whole blood
        dominates; brain tumor is the largest non-blood class.
      </Text>
      <Card>
        <CardHeader>Top tissues</CardHeader>
        <CardBody>
          <BarChart
            categories={TISSUES.map((t) => t.tissue)}
            series={[{ name: "Samples", data: TISSUES.map((t) => t.n) }]}
            height={360}
          />
        </CardBody>
      </Card>
      <Stack gap={8}>
        <H2>Brain region catalogue (Hub)</H2>
        <Table
          headers={["Region", "N"]}
          rows={BRAIN.map((b) => [b.region, b.n.toLocaleString()])}
        />
        <Text tone="secondary">
          No brain_head — region catalogue, not case/control.
        </Text>
      </Stack>
    </Stack>
  );
}

function DiseasesPanel() {
  const top15 = DISEASES.slice(0, 15);
  return (
    <Stack gap={20}>
      <Callout tone="info">
        Disease pack is the Hub nine-pack disease family. Patient fraction =
        cases / pack. Diagnoses come from phenotype_value on case rows;
        controls usually have no diagnosis string. Train only on label_status
        in {"{case, control}"}.
      </Callout>

      <Grid columns={4} gap={12}>
        <Stat value={DISEASE.nPack.toLocaleString()} label="Disease pack" />
        <Stat value={DISEASE.nCase.toLocaleString()} label="Cases (patients)" />
        <Stat
          value={`${(100 * DISEASE.fracOfPack).toFixed(1)}%`}
          label="Patients of pack"
        />
        <Stat
          value={`${(100 * DISEASE.fracOfNinePack).toFixed(1)}%`}
          label="Patients of nine-pack"
        />
      </Grid>

      <Grid columns={3} gap={12}>
        <Stat value={DISEASE.nControl.toLocaleString()} label="Controls" />
        <Stat value={DISEASE.nAdjacent.toLocaleString()} label="Adjacent normal" />
        <Stat value={String(DISEASE.nLabels)} label="Distinct diagnoses" />
      </Grid>

      <Card>
        <CardHeader>Most frequent diagnoses (top 15 cases)</CardHeader>
        <CardBody>
          <BarChart
            horizontal
            categories={top15.map((d) => d.disease)}
            series={[{ name: "Case samples", data: top15.map((d) => d.n) }]}
            height={420}
          />
        </CardBody>
      </Card>

      <Stack gap={8}>
        <H2>All {DISEASE.nLabels} disease labels</H2>
        <Table
          headers={["Disease", "Cases", "% of cases"]}
          rows={DISEASES.map((d) => [
            d.disease,
            d.n.toLocaleString(),
            `${d.pct}%`,
          ])}
        />
      </Stack>

      <Card>
        <CardHeader>Cancer pack (separate Hub family)</CardHeader>
        <CardBody>
          <Grid columns={4} gap={12}>
            <Stat value={CANCER.nPack.toLocaleString()} label="Cancer pack" />
            <Stat value={CANCER.nCase.toLocaleString()} label="Cases" />
            <Stat value={CANCER.nControl.toLocaleString()} label="Controls" />
            <Stat
              value={`${(100 * CANCER.fracOfPack).toFixed(1)}%`}
              label="Patients of pack"
            />
          </Grid>
          <Text tone="secondary">
            Adjacent normal: {CANCER.nAdjacent.toLocaleString()}. Cancer
            diagnoses are a separate census (not listed above).
          </Text>
        </CardBody>
      </Card>
    </Stack>
  );
}
