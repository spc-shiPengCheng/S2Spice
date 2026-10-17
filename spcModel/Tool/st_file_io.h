#ifndef ST_FILE_IO_H
#define ST_FILE_IO_H

#include <QString>

enum class DATAORMODEL
{
    DATA = 0,
    MODEL,
};


void readSparameter(QString filePath,struct snp& data);

void writeSparameter(QString filePath,struct snp& data,DATAORMODEL dataOrModel);

#endif // ST_FILE_IO_H
