#ifndef ST_LOGGER_H
#define ST_LOGGER_H

#include <QObject>

class StLogger: public QObject
{
    Q_OBJECT
public:
    static StLogger* instance();

    void log(const QString& message);

signals:
    void logMessage(const QString& message);

private:
    explicit StLogger(QObject *parent = nullptr);
    static StLogger* m_instance;
};

#endif // ST_LOGGER_H
