# This script is inserted into your .bashrc file when `conda init bash` is run.
# `conda init` is required to be able to use `conda activate [env]`, but then you must open a new shell or reload your .basrc file.
# But Github Actions does not allow you to edit the .bashrc and have changes persist between step shells.
# So we need to run this explicitly to bypass `conda init`.
CONDA_INITIALIZE() {
  echo "CONDA_INITIALIZE [BEGIN]"
  local SHELL=$1
  # echo "$(${CONDA}/bin/conda 'shell.bash' 'hook')"
  __conda_setup="$(${CONDA}/bin/conda shell.${SHELL} 'hook' 2> /dev/null)"
  if [ $? -eq 0 ]; then
    eval "$__conda_setup"
  else 
    if [ -f "${CONDA}/etc/profile.d/conda.sh" ]; then
      . "${CONDA}/etc/profile.d/conda.sh"
    else
      export PATH="${CONDA}/bin:$PATH"
    fi
  fi
  unset __conda_setup
  echo "CONDA_INITIALIZE [END]"
}
