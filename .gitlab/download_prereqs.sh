CHECKOUT_ROOT=$(git rev-parse --show-toplevel)

export XCAP_SPACKAGES_CHECKOUT="${CHECKOUT_ROOT}/extern/xcap_spackages"
export PORTS_OF_CALL_CHECKOUT="${CHECKOUT_ROOT}/extern/ports-of-call"

mkdir -p $(dirname ${XCAP_SPACKAGES_CHECKOUT})
mkdir -p $(dirname ${PORTS_OF_CALL_CHECKOUT})
REPO_URL=$(git remote get-url origin)

if [ ! -d "${XCAP_SPACKAGES_CHECKOUT}" ]; then
  git clone "${REPO_URL%/*/*}/spackages.git" "${XCAP_SPACKAGES_CHECKOUT}"
fi

if [ ! -d "${PORTS_OF_CALL_CHECKOUT}" ]; then
  git clone "${REPO_URL%/*}/ports-of-call.git" "${PORTS_OF_CALL_CHECKOUT}"
fi

if [ -n "${XCAP_SPACKAGES_MR}" ]; then
  git -C ${XCAP_SPACKAGES_CHECKOUT} fetch origin merge-requests/${XCAP_SPACKAGES_MR}/head:mr-${XCAP_SPACKAGES_MR}
  export XCAP_SPACKAGES_REF="${XCAP_SPACKAGES_REF:-"mr-${XCAP_SPACKAGES_MR}"}"
else
  export XCAP_SPACKAGES_REF="${XCAP_SPACKAGES_REF:-main}"
fi

git -C ${XCAP_SPACKAGES_CHECKOUT} checkout ${XCAP_SPACKAGES_REF}

if [ -n "${PORTS_OF_CALL_MR}" ]; then
  git -C ${PORTS_OF_CALL_CHECKOUT} fetch origin merge-requests/${PORTS_OF_CALL_MR}/head:mr-${PORTS_OF_CALL_MR}
  export PORTS_OF_CALL_REF="${PORTS_OF_CALL_REF:-"mr-${PORTS_OF_CALL_MR}"}"
else
  export PORTS_OF_CALL_REF="${PORTS_OF_CALL_REF:-main}"
fi

git -C ${PORTS_OF_CALL_CHECKOUT} checkout ${PORTS_OF_CALL_REF}
