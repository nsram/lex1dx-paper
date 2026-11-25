# lex_strict_fp(<target>): compile <target> with IEEE-754 semantics and no floating-point
# contraction (no FMA fusion of a*b+c), so results are bit-identical on every platform and the
# tie rule behaves the same everywhere.
#
#   GCC, Clang, AppleClang : -ffp-contract=off (Apple clang and GCC otherwise fuse on arm64 / FMA hardware)
#   MSVC (VS 2022, 19.30+) : /fp:precise, which does not generate contractions from VS 2022 on
#                            (contractions are opt-in there via /fp:contract; there is no /fp:contract-)
#   MSVC before VS 2022    : /fp:strict, because /fp:precise could contract in those versions
#   clang-cl               : /fp:precise plus /clang:-ffp-contract=off
# See https://learn.microsoft.com/cpp/build/reference/fp-specify-floating-point-behavior
function(lex_strict_fp target)
  if(MSVC)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      target_compile_options(${target} PRIVATE /fp:precise /clang:-ffp-contract=off)
    elseif(MSVC_VERSION LESS 1930)
      target_compile_options(${target} PRIVATE /fp:strict)
    else()
      target_compile_options(${target} PRIVATE /fp:precise)
    endif()
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${target} PRIVATE -ffp-contract=off)
  else()
    message(WARNING "lex: unknown compiler ${CMAKE_CXX_COMPILER_ID}; make sure it does not contract a*b+c into FMA")
  endif()
endfunction()
