// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef __HERMES_DATA_TABLE_H
#define __HERMES_DATA_TABLE_H

// data table row
struct DataTableRow
{
    double key;
    double value;
    DataTableRow *next;
};

class DataTable
{
public:
    DataTable();
    ~DataTable();

    void clear();
    void remove(double key);

    void add(double key, double value);
    void add(double *key, double *value, int count);

    int size();
    double min_key();
    double max_key();
    double min_value();
    double max_value();

    double value(double key);
    double derivative(double key);

    void print();
    void save(const char *filename, double start, double end, int count);

private:
    DataTableRow *m_data;
};

#endif
