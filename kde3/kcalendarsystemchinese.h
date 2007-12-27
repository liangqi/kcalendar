/*
    Copyright (c) 2007 Liang Qi <cavendish.qi@gmail.com>
        Calendar conversion routines based on:
        1. ccal v2.4 by Zhuo Meng <zhuo@thunder.cwru.edu>
        2. NOVAS-C v2.0 (1 Nov 98) by
               U. S. Naval Observatory
               Astronomical Applications Dept.
               3450 Massachusetts Ave., NW
               Washington, DC  20392-5420
        3. Lunar Outreach Services by Christopher Osburn(1996)
 
    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.
 
    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Library General Public License for more details.
 
    You should have received a copy of the GNU Library General Public License
    along with this library; see the file COPYING.LIB.  If not, write to
    the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
    Boston, MA 02110-1301, USA.
*/

#ifndef KCALENDARSYSTEMCHINESE_H
#define KCALENDARSYSTEMCHINESE_H

#include <qdatetime.h>
#include <qstring.h>

#include "kcalendarsystem.h"

class KCalendarSystemChinesePrivate;

/**
 * @internal
 * This is the Chinese calendar implementation.
 *
 * The Chinese calendar is the traditional calendar used by Chinese
 *
 * @see KLocale,KCalendarSystem,KCalendarSystemFactory
 *
 * @author Liang Qi <cavendish.qi@gmail.com>
 * @since 3.5.9
 */
class KDECORE_EXPORT KCalendarSystemChinese : public KCalendarSystem
{
public:
  /** Constructor. Just like KCalendarSystem::KCalendarSystem(). */
  explicit KCalendarSystemChinese( const KLocale * locale = 0 );
  virtual ~KCalendarSystemChinese();

  virtual int year (const QDate & date) const;
  virtual int month (const QDate & date) const;
  virtual int day (const QDate & date) const;
  virtual int dayOfWeek (const QDate & date) const;
  virtual int dayOfYear (const QDate & date) const;

  virtual bool setYMD(QDate & date, int y, int m, int d) const;

  virtual QDate addYears(const QDate & date, int nyears) const;
  virtual QDate addMonths(const QDate & date, int nmonths) const;
  virtual QDate addDays(const QDate & date, int ndays) const;

  virtual int monthsInYear (const QDate & date) const;
  virtual int daysInYear (const QDate & date) const;
  virtual int daysInMonth (const QDate & date) const;
  virtual int weeksInYear(int year) const;
  virtual int weekNumber(const QDate& date, int * yearNum = 0) const;

  virtual QString monthName (int month, int year, bool shortName = false) const;
  virtual QString monthName (const QDate & date, bool shortName = false ) const;
  virtual QString monthNamePossessive(int month, int year, bool shortName = false) const;
  virtual QString monthNamePossessive(const QDate & date, bool shortName = false ) const;
  virtual QString weekDayName (int weekDay, bool shortName = false) const;
  virtual QString weekDayName (const QDate & date, bool shortName = false) const;

  virtual QString dayString(const QDate & pDate, bool bShort) const;
  virtual QString yearString(const QDate & pDate, bool bShort) const;
  virtual int dayStringToInteger(const QString & sNum, int & iLength) const;
  virtual int yearStringToInteger(const QString & sNum, int & iLength) const;

  virtual int minValidYear () const;
  virtual int maxValidYear () const;
  virtual int weekDayOfPray () const;

  virtual QString calendarName() const;

  virtual bool isLunar() const;
  virtual bool isLunisolar() const;
  virtual bool isSolar() const;

private:
    KCalendarSystemChinesePrivate * d;
};

#endif // KCALENDARSYSTEMCHINESE_H
