#ifndef __ACCENV_PROF_ACC_H
#define __ACCENV_PROF_ACC_H

void prof_acc_start_capture();
void prof_acc_stop_capture();
void prof_acc_begin(const char* name);
void prof_acc_end(const char* name);

#endif

