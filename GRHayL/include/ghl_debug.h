#ifndef GHL_DEBUG_H_
#define GHL_DEBUG_H_

#define _indent_1(string) (int)((n_subdiv - strlen(string) - 1)/2)
#define _indent_2(string) (int)(n_subdiv - 1 - _indent_1(string))
#define _print_6(a, b, c, d, e, f)                                       \
  fprintf(stderr, "|%*s%-*s|%*s%-*s|%*s%-*s|%*s%-*s|%*s%-*s|%*s%-*s|\n", \
          _indent_1(a), "", _indent_2(a), a,                             \
          _indent_1(b), "", _indent_2(b), b,                             \
          _indent_1(c), "", _indent_2(c), c,                             \
          _indent_1(d), "", _indent_2(d), d,                             \
          _indent_1(e), "", _indent_2(e), e,                             \
          _indent_1(f), "", _indent_2(f), f);

#define _print_7(a, b, c, d, e, f, g)                                            \
  fprintf(stderr, "|%*s%-*s|%*s%-*s|%*s%-*s|%*s%-*s|%*s%-*s|%*s%-*s|%*s%-*s|\n", \
          _indent_1(a), "", _indent_2(a), a,                                     \
          _indent_1(b), "", _indent_2(b), b,                                     \
          _indent_1(c), "", _indent_2(c), c,                                     \
          _indent_1(d), "", _indent_2(d), d,                                     \
          _indent_1(e), "", _indent_2(e), e,                                     \
          _indent_1(f), "", _indent_2(f), f,                                     \
          _indent_1(g), "", _indent_2(g), g);

static inline void ghl_debug_print_prims(
      const ghl_primitive_quantities *restrict p) {

  fprintf(stderr,
".------------------------.------------------------.------------------------.------------------------.------------------------.------------------------.\n"
"|          rho           |        Pressure        |          eps           |        Entropy         |          Y_e           |      Temperature       |\n"
"| %22.15e | %22.15e | %22.15e | %22.15e | %22.15e | %22.15e |\n"
".------------------------.------------------------.------------------------.------------------------.------------------------.------------------------.\n"
"|          VU0           |          VU1           |          VU2           |          BU0           |          BU1           |          BU2           |\n"
"| %22.15e | %22.15e | %22.15e | %22.15e | %22.15e | %22.15e |\n"
".------------------------.------------------------.------------------------.------------------------.------------------------.------------------------.\n\n",
  p->rho, p->press, p->eps, p->entropy, p->Y_e, p->temperature,
  p->vU[0], p->vU[1], p->vU[2], p->BU[0], p->BU[1], p->BU[2]);
}

static inline void ghl_debug_print_cons(
      const ghl_conservative_quantities *restrict c) {

  fprintf(stderr,
".------------------------.------------------------.------------------------.------------------------.------------------------.------------------------.------------------------.\n"
"|           D            |          tau           |          SD0           |          SD1           |          SD2           |        Entropy         |          Y_e           |\n"
"| %22.15e | %22.15e | %22.15e | %22.15e | %22.15e | %22.15e | %22.15e |\n"
".------------------------.------------------------.------------------------.------------------------.------------------------.------------------------.------------------------.\n\n",
          c->tau, c->SD[0], c->SD[1], c->SD[2], c->rho, c->entropy, c->Y_e);
}

#define ghl_debug_print_primitives(p) {                                \
  fprintf(stderr, "DEBUG (Primitive Quantities) func: %s, line: %d\n", \
          __func__, __LINE__);                                         \
  ghl_debug_print_prims(p);                                            \
}

#define ghl_debug_print_conservatives(c) {                                \
  fprintf(stderr, "DEBUG (Conservative Quantities) func: %s, line: %d\n", \
          __func__, __LINE__);                                            \
  ghl_debug_print_cons(c);                                                \
}

#endif // GHL_DEBUG_H_
