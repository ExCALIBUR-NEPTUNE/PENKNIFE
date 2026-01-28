#!/usr/bin/env bash

_run_eg_completions() {
    local REPO_ROOT=$( cd -- "$(realpath $( dirname -- "${BASH_SOURCE[0]}" )/..)" &> /dev/null && pwd )
    local BASE_DIR=$REPO_ROOT/examples

    local cur prev
    local level1 level2 level3

    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    level1="${COMP_WORDS[1]}"
    level2="${COMP_WORDS[2]}"
    level3="${COMP_WORDS[3]}"

    case "${COMP_CWORD}" in
        1)
            # Complete top-level directories
            COMPREPLY=($(compgen -W "$(ls -1d "${BASE_DIR}"/*/ 2>/dev/null | xargs -n1 basename)" -- "$cur"))
            ;;
        2)
            # Complete directories inside level1
            if [[ -d "${BASE_DIR}/${level1}" ]]; then
                COMPREPLY=($(compgen -W "$(ls -1d "${BASE_DIR}/${level1}"/*/ 2>/dev/null | xargs -n1 basename)" -- "$cur"))
            fi
            ;;
        3)
            # Complete directories inside level2
            if [[ -d "${BASE_DIR}/${level1}/${level2}" ]]; then
                COMPREPLY=($(compgen -W "$(ls -1d "${BASE_DIR}/${level1}/${level2}"/*/ 2>/dev/null | xargs -n1 basename)" -- "$cur"))
            fi
            ;;
        4)
            # Complete directories inside level3
            if [[ -d "${BASE_DIR}/${level1}/${level2}/${level3}" ]]; then
                COMPREPLY=($(compgen -W "$(ls -1d "${BASE_DIR}/${level1}/${level2}/${level3}"/*/ 2>/dev/null | xargs -n1 basename)" -- "$cur"))
            fi
            ;;
    esac
}

complete -F _run_eg_completions ./scripts/run_eg.sh
