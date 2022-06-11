#ifndef ANALYZE_H
#define ANALYZE_H
double rgsq(PARAMS * p, double * coords, int t);
int ** unique_paths(PARAMS * p, double ** coords_arr);
double min_rgsq(PARAMS * p, double * coords);
#endif// ANALYZE_H
