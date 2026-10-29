#!/usr/bin/env bash
set -euo pipefail

# Par defaut : nombre de coeurs disponibles (ou 4 si nproc non present)
MAX_JOBS="$(nproc 2>/dev/null || echo 4)"

# Gestion des arguments optionnels (-j N) et du mode
while [[ $# -gt 0 ]]; do
  case "$1" in
    -j|--jobs)
      MAX_JOBS="$2"
      shift 2
      ;;
    pilot|reference|full)
      mode="$1"
      shift
      ;;
    *)
      echo "Usage : $0 [-j NUM_JOBS] {pilot|reference|full}" >&2
      exit 2
      ;;
  esac
done

if [[ -z "${mode:-}" ]]; then
  echo "Usage : $0 [-j NUM_JOBS] {pilot|reference|full}" >&2
  exit 2
fi

echo "Mode: $mode | Concurrence max: $MAX_JOBS"

project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$project_dir"
qcc_bin="${QCC:-$(command -v qcc || true)}"
if [[ -z "$qcc_bin" && -x "${BASILISK:-$HOME/basilisk/src}/qcc" ]]; then
  qcc_bin="${BASILISK:-$HOME/basilisk/src}/qcc"
fi
if [[ -z "$qcc_bin" ]]; then
  echo "qcc introuvable : ajouter Basilisk au PATH ou definir QCC." >&2
  exit 1
fi
for tool in gcc python3 setsid sha256sum; do
  if ! command -v "$tool" >/dev/null; then
    echo "$tool introuvable." >&2
    exit 1
  fi
done
python3 -c 'import numpy, matplotlib'
if ! command -v jupyter >/dev/null; then
  echo "Jupyter absent : necessaire seulement pour ouvrir analyse_campagne.ipynb."
fi
export OMP_NUM_THREADS=1
source_hash="$(sha256sum jurin_axi_amr.c | cut -d ' ' -f 1)"

# Tableaux pour gerer les processus actifs en parallele
declare -A ACTIVE_PIDS=()
declare -A ACTIVE_DIRS=()

finish () {
  local rc=$?
  trap - EXIT INT TERM
  # Arreter tous les groupes de processus enfants encore actifs
  for pid in "${!ACTIVE_PIDS[@]}"; do
    kill -TERM -- "-$pid" 2>/dev/null || true
  done
  wait 2>/dev/null || true
  for dir in "${ACTIVE_DIRS[@]}"; do
    if [[ -d "$dir" ]]; then
      if (( rc == 130 || rc == 143 )); then
        echo interrupted > "$dir/status"
      else
        echo failed > "$dir/status"
      fi
    fi
  done
  exit "$rc"
}
trap finish EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

validate_result () {
  python3 - "$1" "$2" "$3" <<'PY'
import json
import sys
from pathlib import Path
import numpy as np

directory, level, tolerance = Path(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
p = json.loads((directory / 'parametres.json').read_text())
expected = dict(minlevel=5, maxlevel=level, f_tol=1e-3, u_tol=tolerance,
                t_end=0.8, t_compare=0.15, length_m=0.08, radius_m=0.005,
                sample_nx=1024, sample_ny=64)
if p != expected:
    raise ValueError(f'Parametres compiles inattendus : {p}')
h = np.loadtxt(directory / 'hauteur_axi_amr.dat', ndmin=2)
c = np.loadtxt(directory / 'cellules.dat', ndmin=2)
v = np.loadtxt(directory / 'vitesse_t015.dat', ndmin=2)
for name, data in [('hauteur', h), ('cellules', c), ('vitesse', v)]:
    if not np.isfinite(data).all() or np.any(np.abs(data) >= 1e29):
        raise ValueError(f'{name} : valeur invalide ou nodata')
if h.shape != (1601, 3) or not np.allclose(h[:, 0], np.arange(1601)*0.0005, rtol=0, atol=1e-10):
    raise ValueError('Les hauteurs ne couvrent pas tout [0, 0.8] s')
if np.any(h[:, 1:] < 0) or np.any(h[:, 1:] > p['length_m']):
    raise ValueError('Hauteur hors du tube')
if c.shape != (81, 2) or not np.allclose(c[:, 0], np.arange(81)*0.01, rtol=0, atol=1e-10):
    raise ValueError('Le comptage ne couvre pas tout [0, 0.8] s')
if (c[0, 1] != 64 or np.any(c[:, 1] < 64)
        or np.any(c[:, 1] > 2**(2*level - 4)) or np.any(c[:, 1] != np.floor(c[:, 1]))):
    raise ValueError('Nombre de cellules physiques incoherent')
nx, ny = p['sample_nx'], p['sample_ny']
xx = np.repeat((np.arange(nx) + 0.5)*p['length_m']/nx, ny)
yy = np.tile((np.arange(ny) + 0.5)*p['radius_m']/ny, nx)
if v.shape != (nx*ny, 3) or not np.array_equal(v[:, :2], np.column_stack((xx, yy))):
    raise ValueError('Grille de vitesse differente de la grille commune')
with (directory / 'vitesse_t015.dat').open() as stream:
    sample_time = float(stream.readline().split('=', 1)[1])
if abs(sample_time - p['t_compare']) > 1e-12:
    raise ValueError('Instant du champ de vitesse incorrect')
timing = dict(line.split() for line in (directory / 'temps.dat').read_text().splitlines())
if not np.isfinite(float(timing['real_s'])) or float(timing['real_s']) < 0:
    raise ValueError('Temps reel invalide')
print('OK : duree complete, cellules actives, grille commune 1024 x 64, valeurs finies.')
PY
}

run_case () {
  local level="$1" tolerance="$2" name="$3"
  local directory="$project_dir/campaign/$name"
  mkdir -p "$directory"
  if [[ -f "$directory/status" && "$(cat "$directory/status")" == complete
        && -f "$directory/source.sha256" && "$(cat "$directory/source.sha256")" == "$source_hash" ]] \
     && validate_result "$directory" "$level" "$tolerance" > "$directory/validation.log" 2>&1; then
    echo "[$name] deja termine et verifie."
    return 0
  fi

  echo running > "$directory/status"
  rm -f -- "$directory"/{hauteur_axi_amr.dat,cellules.dat,vitesse_t015.dat,parametres.json,temps.dat,validation.log,jurin}
  echo "$source_hash" > "$directory/source.sha256"
  echo "[$name] compilation (MAXLEVEL=$level, U_TOL=$tolerance m/s)..."
  local -a build=("$qcc_bin" -autolink -disable-dimensions -O2 -Wall
    -DMINLEVEL=5 "-DMAXLEVEL=$level" -DF_TOL=1e-3 "-DU_TOL=$tolerance"
    -DJURIN_EXPORT_VIDEO=0 jurin_axi_amr.c -o "$directory/jurin" -lm)
  printf '%q ' "${build[@]}" > "$directory/compile.log"
  printf '\n' >> "$directory/compile.log"
  if ! "${build[@]}" >> "$directory/compile.log" 2>&1; then
    echo failed > "$directory/status"
    echo "[$name] echec de compilation (voir $directory/compile.log)" >&2
    return 1
  fi

  echo "[$name] calcul en cours jusqu'a 0.8 s..."
  local started
  started="$(python3 -c 'import time; print(time.monotonic())')"

  # Lance jurin dans un groupe de processus dedie
  if [[ -x /usr/bin/time ]]; then
    (cd "$directory"; exec setsid /usr/bin/time -f 'real_s %e\nmax_rss_kb %M' -o temps.dat ./jurin) \
      > "$directory/simulation.log" 2>&1
  else
    (cd "$directory"; exec setsid ./jurin) > "$directory/simulation.log" 2>&1
  fi

  if [[ ! -x /usr/bin/time ]]; then
    python3 - "$started" > "$directory/temps.dat" <<'PY'
import sys
import time
print(f'real_s {time.monotonic() - float(sys.argv[1]):.6f}')
print('max_rss_kb nan')
PY
  fi

  if ! validate_result "$directory" "$level" "$tolerance" > "$directory/validation.log" 2>&1; then
    echo failed > "$directory/status"
    echo "[$name] echec de validation (voir $directory/validation.log)" >&2
    return 1
  fi

  echo complete > "$directory/status"
  echo "[$name] termine et verifie ($(head -n 1 "$directory/temps.dat"))."
}

# Fonction utilitaire pour vider les taches terminees
cleanup_finished_jobs () {
  for pid in "${!ACTIVE_PIDS[@]}"; do
    if ! kill -0 "$pid" 2>/dev/null; then
      wait "$pid" || true
      unset "ACTIVE_PIDS[$pid]"
      unset "ACTIVE_DIRS[$pid]"
    fi
  done
}

mkdir -p campaign

if [[ "$mode" == reference ]]; then
  run_case 12 1e-5 reference
else
  if [[ "$mode" == pilot ]]; then
    levels=(8 9)
    tolerances=(1e-2 1e-3)
  else
    levels=(6 7 8 9 10 11 12)
    tolerances=(1e-2 3e-3 1e-3 3e-4 1e-4)
  fi

  for level in "${levels[@]}"; do
    for tolerance in "${tolerances[@]}"; do
      printf -v name 'L%02d_U%s' "$level" "$tolerance"

      # Attente active si le nombre maximum de jobs est atteint
      while (( ${#ACTIVE_PIDS[@]} >= MAX_JOBS )); do
        cleanup_finished_jobs
        if (( ${#ACTIVE_PIDS[@]} >= MAX_JOBS )); then
          sleep 0.5
        fi
      done

      # Lancement asynchrone du cas
      run_case "$level" "$tolerance" "$name" &
      sub_pid=$!
      ACTIVE_PIDS["$sub_pid"]=1
      ACTIVE_DIRS["$sub_pid"]="$project_dir/campaign/$name"
    done
  done

  # Attente de la fin de tous les calculs restants
  while (( ${#ACTIVE_PIDS[@]} > 0 )); do
    cleanup_finished_jobs
    sleep 0.5
  done
fi