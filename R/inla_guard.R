# Internal helpers that protect the R session and the machine from a runaway INLA
# process (Unix/macOS only; not exported).
#
# Why: on some platforms the INLA executable can keep running after a numerical
# failure and leak memory at a steady rate (observed: R-INLA 26.08.07, macOS arm64,
# about 30 MB/s, i.e. 100+ GB in an hour). On macOS the leaked memory is compressed,
# so the usual resident-memory figure does not show it. The guard
#   * stops INLA when its real memory use exceeds a limit (physical footprint on
#     macOS, resident memory on Linux),
#   * stops INLA when a model takes longer than the time limit,
#   * stops INLA when the R process that started it has disappeared, and
#   * only ever touches INLA processes that belong to this R process, so parallel
#     jobs of the same user are not affected.
#
# Limits can be changed with options(), e.g. options(MethylCal.inla_max_mem_mb = 16384):
#   MethylCal.inla_timeout    seconds allowed per model       (default 300)
#   MethylCal.inla_max_mem_mb memory allowed per INLA process (default 4096)

.inla_guard_script <- '
R_PID=$1; TIMEOUT=$2; MAXMB=$3; FLAG=$4
OS=$(uname)
# pids of INLA executables that descend from R_PID
inla_children() {
  ps -axo pid=,ppid=,comm= | awk -v root="$R_PID" \'
    { ids[NR] = $1; ppid[$1] = $2; comm[$1] = $3 }
    END { for (i = 1; i <= NR; i++) { x = ids[i]
            if (comm[x] !~ /(^|\\/)inla(\\.mkl)?(\\.run)?$/) continue
            y = x; while (y > 1) { if (ppid[y] == root) { print x; break } y = ppid[y] } } }\'
}
is_inla() { ps -o comm= -p "$1" 2>/dev/null | grep -Eq "(^|/)inla(\\.mkl)?(\\.run)?$"; }
# real memory use of a process in MB
mem_mb() {
  if [ "$OS" = Darwin ]; then
    footprint "$1" 2>/dev/null | awk \'/Footprint:/ { for (i = 1; i <= NF; i++) if ($i == "Footprint:") { v = $(i+1); u = $(i+2) }
         if (u ~ /GB/) v *= 1024; else if (u ~ /KB/) v /= 1024; else if (u ~ /^B/) v /= 1048576; printf "%d", v; exit }\'
  else
    ps -o rss= -p "$1" 2>/dev/null | awk \'{ printf "%d", $1 / 1024 }\'
  fi
}
kill_inla() { for p in $SEEN $(inla_children); do is_inla "$p" && kill -9 "$p" 2>/dev/null; done; }
END=$(( $(date +%s) + TIMEOUT )); SEEN=""; TRIPPED=0
while :; do
  # R is gone: no INLA process may outlive it
  if ! kill -0 "$R_PID" 2>/dev/null; then kill_inla; exit 0; fi
  for p in $(inla_children); do
    case " $SEEN " in *" $p "*) ;; *) SEEN="$SEEN $p" ;; esac
    m=$(mem_mb "$p")
    if [ "${m:-0}" -gt "$MAXMB" ]; then echo memory > "$FLAG"; TRIPPED=1; fi
  done
  if [ "$(date +%s)" -ge "$END" ]; then echo time > "$FLAG"; TRIPPED=1; fi
  # once a limit was hit, also stop the automatic re-run that INLA starts
  if [ "$TRIPPED" = 1 ]; then kill_inla; fi
  sleep 2
done
'

# Start the guard for the current R process. Returns the process and the flag file
# in which the guard records why it stopped INLA ("time" or "memory").
.inla_guard_start <- function(timeout, max_mem_mb) {
  flag   <- tempfile(fileext = ".guard")
  script <- tempfile(fileext = ".sh")
  writeLines(.inla_guard_script, script)
  args   <- c(script, Sys.getpid(), timeout, max_mem_mb, flag)
  perl   <- Sys.which("perl")
  # Start in a new session (setsid) so that the guard survives if the whole
  # process group of the R session is killed; plain sh if perl is not available.
  proc <- if (nzchar(perl))
    processx::process$new(perl, c("-MPOSIX", "-e", "POSIX::setsid(); exec @ARGV", "sh", args))
  else
    processx::process$new("sh", args)
  list(proc = proc, flag = flag, script = script)
}

.inla_guard_stop <- function(guard) {
  if (guard$proc$is_alive()) guard$proc$kill()
  unlink(c(guard$flag, guard$script), force = TRUE)
  invisible(NULL)
}

# Why did the guard stop INLA? "time", "memory" or "" (it did not).
.inla_guard_reason <- function(guard) {
  if (!file.exists(guard$flag)) return("")
  trimws(paste(readLines(guard$flag, warn = FALSE), collapse = ""))
}

# Kill INLA processes started by this R process (and only those).
.inla_kill_children <- function() {
  ps <- suppressWarnings(system("ps -axo pid=,ppid=,comm=", intern = TRUE, ignore.stderr = TRUE))
  if (!length(ps)) return(invisible(NULL))
  m   <- regmatches(ps, regexec("^\\s*([0-9]+)\\s+([0-9]+)\\s+(\\S+)", ps))
  m   <- m[lengths(m) == 4]
  pid <- as.integer(vapply(m, `[`, "", 2)); ppid <- as.integer(vapply(m, `[`, "", 3)); comm <- vapply(m, `[`, "", 4)
  up  <- stats::setNames(ppid, pid)
  inla <- pid[grepl("(^|/)inla(\\.mkl)?(\\.run)?$", comm)]
  me  <- Sys.getpid()
  for (p in inla) {
    q <- p
    while (!is.na(q) && q > 1) { if (identical(unname(up[as.character(q)]), me)) { tools::pskill(p, tools::SIGKILL); break }
                                  q <- unname(up[as.character(q)]) }
  }
  invisible(NULL)
}
