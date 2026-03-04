import re
import sys

# Regex for lines like:
# else if 0100<=&SIC_Code<=0199 then do; FFI48=1; FFI48_desc='Agric'; end;
pat_range = re.compile(
    r"else\s+if\s+(\d+)\s*<=\s*&\w+\s*<=\s*(\d+)\s*then\s+do;\s*"
    r"FFI48\s*=\s*(\d+)\s*;\s*FFI48_desc\s*=\s*'([^']+)'\s*;\s*end;",
    re.IGNORECASE
)

# Regex for "if missing(&SIC_Code) then FFI48=.;" (optional)
pat_missing = re.compile(r"if\s+missing\(\s*&\w+\s*\)\s+then\s+FFI48\s*=\s*\.\s*;", re.IGNORECASE)

def main(path):
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    out = []
    out.append("* ---- BEGIN GENERATED RULES ----")
    out.append("* Assumes tempvar sicnum exists and output vars `gen' and `desc' exist.")
    out.append("* Guarded with missing(`gen') to mimic SAS else-if (first match wins).")
    out.append("")

    for line in lines:
        line = line.strip()

        # Skip comment blocks and empty lines
        if not line or line.startswith("/*") or line.startswith("*"):
            continue

        if pat_missing.search(line):
            # In Stata we typically just initialize to missing, so nothing required.
            continue

        m = pat_range.search(line)
        if not m:
            continue

        lo, hi, code, desc = m.group(1), m.group(2), m.group(3), m.group(4)

        # Convert to ints to drop leading zeros safely
        lo_i, hi_i, code_i = int(lo), int(hi), int(code)

        if lo_i == hi_i:
            cond = f"`sicnum'=={lo_i}"
        else:
            cond = f"inrange(`sicnum', {lo_i}, {hi_i})"

        out.append(f"replace `gen'  = {code_i}     if `touse' & missing(`gen') & {cond}")
        out.append(f"replace `desc' = \"{desc}\" if `touse' & `gen'=={code_i} & `desc'==\"\"")
        out.append("")

    out.append("* ---- END GENERATED RULES ----")

    print("\n".join(out))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python sas_ffi48_to_stata.py path/to/ffi48.sas", file=sys.stderr)
        sys.exit(2)
    main(sys.argv[1])