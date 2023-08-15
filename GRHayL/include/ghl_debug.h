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

static inline
void ghl_debug_print_prims( const ghl_primitive_quantities *restrict p ) {
  const unsigned short n_subdiv = 25;
  const unsigned short n_div    = 6;
  const unsigned short n_total  = n_subdiv * n_div;
  char subdiv[n_subdiv], div[n_total+3];

  // Prepare sub-division
  subdiv[0] = '.';
  for(int i=1;i<n_subdiv;i++) subdiv[i] = '-';

  // Prepare division
  div[0] = '\0';
  for(int i=0;i<n_div;i++) sprintf(div, "%s%s", div, subdiv);
  div[n_total  ] = '.';
  div[n_total+1] = '\n';
  div[n_total+2] = '\0';

  fputs(div, stderr);
  _print_6("VU0", "VU1", "VU2", "BU0", "BU1", "BU2");
  fputs(div, stderr);
  fprintf(stderr, "| %22.15e | %22.15e | %22.15e | %22.15e | %22.15e | %22.15e |\n",
          p->vU[0], p->vU[1], p->vU[2], p->BU[0], p->BU[1], p->BU[2]);
  fputs(div, stderr);
  _print_6("rho", "Pressure", "eps", "Entropy", "Y_e", "Temperature");
  fputs(div, stderr);
  fprintf(stderr, "| %22.15e | %22.15e | %22.15e | %22.15e | %22.15e | %22.15e |\n",
          p->rho, p->press, p->eps, p->entropy, p->Y_e, p->temperature);
  fputs(div, stderr);
}

static inline
void ghl_debug_print_cons( const ghl_conservative_quantities *restrict c ) {
  const unsigned short n_subdiv = 25;
  const unsigned short n_div    = 7;
  const unsigned short n_total  = n_subdiv * n_div;
  char subdiv[n_subdiv], div[n_total+3];

  // Prepare sub-division
  subdiv[0] = '.';
  for(int i=1;i<n_subdiv;i++) subdiv[i] = '-';

  // Prepare division
  div[0] = '\0';
  for(int i=0;i<n_div;i++) sprintf(div, "%s%s", div, subdiv);
  div[n_total  ] = '.';
  div[n_total+1] = '\n';
  div[n_total+2] = '\0';

  fputs(div, stderr);
  _print_7("tau", "SD0", "SD1", "SD2", "rho", "Entropy", "Y_e");
  fputs(div, stderr);
  fprintf(stderr, "| %22.15e | %22.15e | %22.15e | %22.15e | %22.15e | %22.15e | %22.15e |\n",
          c->tau, c->SD[0], c->SD[1], c->SD[2], c->rho, c->entropy, c->Y_e);
  fputs(div, stderr);
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
