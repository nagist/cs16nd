//===== Copyright ?1996-2005, Valve Corporation, All rights reserved. ======//
//
// Purpose: 
//
// $Workfile:     $
// $Date:         $
//
//-----------------------------------------------------------------------------
// $Log: $
//
// $NoKeywords: $
//===========================================================================//

class CPerformanceCounter
{
public:
	void InitializePerformanceCounter(void);
	float GetCurTime(void);

public:
	int m_iLowShift;
	float m_flPerfCounterFreq;
	float m_flCurrentTime;
	float m_flLastCurrentTime;
};