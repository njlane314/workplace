#!/usr/bin/env bash

_heron()
{
  local cur prev
  COMPREPLY=()
  cur="${COMP_WORDS[COMP_CWORD]}"
  prev="${COMP_WORDS[COMP_CWORD-1]}"

  local commands="art sample event macro paths env help -h --help"


  _heron_has_macros()
  {
    local root="$1"
    [[ -d "${root}/macros/plot/macro" || -d "${root}/macros/evd/macro" || -d "${root}/macros/standalone/macro" || -d "${root}/macros/io/macro" || -d "${root}/macro_packs/default/plot/macro" || -d "${root}/macro_packs/default/evd/macro" || -d "${root}/macro_packs/default/standalone/macro" || -d "${root}/macro_packs/default/io/macro" ]]
  }

  _heron_macro_root()
  {
    local root="$1"
    if [[ -d "${root}/macros/plot/macro" || -d "${root}/macros/evd/macro" || -d "${root}/macros/standalone/macro" || -d "${root}/macros/io/macro" ]]; then
      printf "%s" "${root}/macros"
      return 0
    fi
    if [[ -d "${root}/macro_packs/default/plot/macro" || -d "${root}/macro_packs/default/evd/macro" || -d "${root}/macro_packs/default/standalone/macro" || -d "${root}/macro_packs/default/io/macro" ]]; then
      printf "%s" "${root}/macro_packs/default"
      return 0
    fi
    return 1
  }

  _heron_find_root()
  {
    local dir
    local base
    local parent
    local exe
    local exe_dir

    if [[ -n "${HERON_REPO_ROOT:-}" ]] && _heron_has_macros "${HERON_REPO_ROOT}"; then
      printf "%s" "${HERON_REPO_ROOT}"
      return 0
    fi

    if [[ -n "${HERON_ROOT:-}" ]] && _heron_has_macros "${HERON_ROOT}"; then
      printf "%s" "${HERON_ROOT}"
      return 0
    fi

    exe="$(command -v heron 2>/dev/null || true)"
    if [[ -n "${exe}" ]]; then
      exe_dir="$(dirname "$(readlink -f "${exe}" 2>/dev/null || printf "%s" "${exe}")")"
      dir="${exe_dir}"
      while [[ -n "${dir}" && "${dir}" != "/" ]]; do
        if _heron_has_macros "${dir}"; then
          printf "%s" "${dir}"
          return 0
        fi
        dir="$(dirname "${dir}")"
      done
    fi

    dir="${PWD}"
    while [[ "${dir}" != "/" ]]; do
      if _heron_has_macros "${dir}"; then
        printf "%s" "${dir}"
        return 0
      fi
      base="$(basename "${dir}")"
      dir="$(dirname "${dir}")"
    done
    return 1
  }

  _heron_list_macros()
  {
    local repo_root
    local macro_dir
    local macro
    local evd_dir

    repo_root="$(_heron_find_root 2>/dev/null || true)"
    if [[ -n "${repo_root}" ]]; then
      local macro_root
      macro_root="$(_heron_macro_root "${repo_root}" 2>/dev/null || true)"
      if [[ -z "${macro_root}" ]]; then
        return 0
      fi

      macro_dir="${macro_root}/plot/macro"
      if [[ -d "${macro_dir}" ]]; then
        for macro in "${macro_dir}"/*.C; do
          if [[ -f "${macro}" ]]; then
            basename "${macro}"
          fi
        done
      fi
      macro_dir="${macro_root}/standalone/macro"
      if [[ -d "${macro_dir}" ]]; then
        for macro in "${macro_dir}"/*.C; do
          if [[ -f "${macro}" ]]; then
            basename "${macro}"
          fi
        done
      fi
      evd_dir="${macro_root}/evd/macro"
      if [[ -d "${evd_dir}" ]]; then
        for macro in "${evd_dir}"/*.C; do
          if [[ -f "${macro}" ]]; then
            basename "${macro}"
          fi
        done
      fi
      macro_dir="${macro_root}/io/macro"
      if [[ -d "${macro_dir}" ]]; then
        for macro in "${macro_dir}"/*.C; do
          if [[ -f "${macro}" ]]; then
            basename "${macro}"
          fi
        done
      fi
      return 0
    fi

    heron macro list 2>/dev/null | awk '/\.C$/ {print $1}'
  }

  if [[ ${COMP_CWORD} -le 1 ]]; then
    COMPREPLY=( $(compgen -W "${commands}" -- "${cur}") )
    return 0
  fi

  if [[ "${prev}" == "help" || "${prev}" == "-h" || "${prev}" == "--help" ]]; then
    COMPREPLY=( $(compgen -W "${commands}" -- "${cur}") )
    return 0
  fi

  if [[ "${COMP_WORDS[1]}" == "macro" ]]; then
    local macros
    macros="$(_heron_list_macros)"
    if [[ ${COMP_CWORD} -eq 2 ]]; then
      COMPREPLY=( $(compgen -W "${macros}" -- "${cur}") )
      return 0
    fi
  fi

  COMPREPLY=()
  return 0
}

complete -F _heron heron
