import fs from "node:fs";
import path from "node:path";

const DEFAULT_INPUT = "/home/mayuan/code/my_script/test/acnn/cam";
const DEFAULT_OUTPUT_DIR = "/home/mayuan/code/my_script/test/acnn";

const inputPath = process.argv[2] ?? DEFAULT_INPUT;
const outputDir = process.argv[3] ?? DEFAULT_OUTPUT_DIR;

const stablePath = path.join(outputDir, "stable.csv");
const unstablePath = path.join(outputDir, "unstable.csv");

function looksLikeFormula(token: string): boolean {
  return /^[A-Z][a-z]?\d*(?:[A-Z][a-z]?\d*)*$/.test(token);
}

function parseFormulaCounts(formula: string): Record<string, number> {
  const counts: Record<string, number> = {};
  const regex = /([A-Z][a-z]?)(\d*)/g;
  let match: RegExpExecArray | null;
  while ((match = regex.exec(formula)) !== null) {
    const element = match[1];
    const numberText = match[2];
    const count = numberText === "" ? 1 : Number(numberText);
    counts[element] = (counts[element] ?? 0) + count;
  }
  return counts;
}

function findFormula(tokens: string[]): string | null {
  for (const token of tokens) {
    if (looksLikeFormula(token)) {
      return token;
    }
  }
  return null;
}

function normalizeEnthalpy(value: number, raw: string): string {
  if (Number.isNaN(value)) {
    return raw;
  }
  if (Math.abs(value) <= 1e-8) {
    return "0";
  }
  return String(value);
}

if (!fs.existsSync(inputPath)) {
  throw new Error(`输入文件不存在: ${inputPath}`);
}

const content = fs.readFileSync(inputPath, "utf8");
const lines = content.split(/\r?\n/);

const header = "Number,formula,H,Lu,Li,enthalpy";
const stableRows: string[] = [header];
const unstableRows: string[] = [header];

const warnings: string[] = [];

for (let index = 0; index < lines.length; index += 1) {
  const rawLine = lines[index];
  if (rawLine.trim() === "") {
    continue;
  }

  const tokens = rawLine.trim().split(/\s+/);
  if (tokens.length < 7) {
    warnings.push(`第 ${index + 1} 行列数不足: ${rawLine}`);
    continue;
  }

  const enthalpyToken = tokens[5];
  const enthalpyValue = Number(enthalpyToken);
  if (Number.isNaN(enthalpyValue)) {
    warnings.push(`第 ${index + 1} 行第 6 列非数字: ${enthalpyToken}`);
    continue;
  }

  const formula = findFormula(tokens.slice(6));
  if (formula === null) {
    warnings.push(`第 ${index + 1} 行未找到化学式: ${rawLine}`);
    continue;
  }

  const counts = parseFormulaCounts(formula);
  const hCount = counts.H ?? 0;
  const luCount = counts.Lu ?? 0;
  const liCount = counts.Li ?? 0;

  const number = index + 1;
  const enthalpyText = normalizeEnthalpy(enthalpyValue, enthalpyToken);
  const row = `${number},${formula},${hCount},${luCount},${liCount},${enthalpyText}`;

  if (Math.abs(enthalpyValue) <= 1e-8) {
    stableRows.push(row);
  } else if (enthalpyValue > 0 && enthalpyValue <= 0.05) {
    unstableRows.push(row);
  }
}

fs.writeFileSync(stablePath, `${stableRows.join("\n")}\n`, "utf8");
fs.writeFileSync(unstablePath, `${unstableRows.join("\n")}\n`, "utf8");

if (warnings.length > 0) {
  const warningText = warnings.join("\n");
  process.stderr.write(`${warningText}\n`);
}
