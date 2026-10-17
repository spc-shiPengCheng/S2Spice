#include "st_logger.h"

StLogger* StLogger::m_instance = nullptr;

StLogger::StLogger(QObject *parent) : QObject(parent)
{
}

StLogger* StLogger::instance()
{
    if (!m_instance) {
        m_instance = new StLogger();
    }
    return m_instance;
}

void StLogger::log(const QString& message)
{
    emit logMessage(message);
}
