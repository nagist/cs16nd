/***
*
*	Copyright (c) 1996-2002, Valve LLC. All rights reserved.
*
*	This product contains software technology licensed from Id
*	Software, Inc. ("Id Technology").  Id Technology (c) 1996 Id Software, Inc.
*	All Rights Reserved.
*
*   Use, distribution, and modification of this source code and/or resulting
*   object code is restricted to non-commercial enhancements to products from
*   Valve LLC.  All other use, distribution, or modification is prohibited
*   without written permission from Valve LLC.
*
****/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#if defined _MSC_VER && _MSC_VER >= 1400
#ifndef _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_DEPRECATE
#endif
#pragma warning(disable:4996)
#endif

#include "../common/nowin.h"

typedef int BOOL;
#define TRUE 1
#define FALSE 0

typedef int qboolean;
typedef unsigned char byte;

#include "../utils/common/mathlib.h"
#include "const.h"
#include "progdefs.h"
#include "edict.h"
#include "eiface.h"
#include "studio.h"

#ifndef ACTIVITY_H
#include "activity.h"
#endif

#include "activitymap.h"

#ifndef ANIMATION_H
#include "animation.h"
#endif

#ifndef SCRIPTEVENT_H
#include "scriptevent.h"
#endif

#ifndef ENGINECALLBACK_H
#include "enginecallback.h"
#endif

#include "r_studioint.h"
#include "com_model.h"

#define EXPORT _declspec(dllexport)

extern globalvars_t *gpGlobals;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ANIM_WALK_SEQUENCE 3
#define ANIM_JUMP_SEQUENCE 6
#define ANIM_SWIM_1 8
#define ANIM_SWIM_2 9
#define ANIM_FIRST_DEATH_SEQUENCE 101

#ifdef _MSC_VER
#pragma warning(disable:4244)
#endif

int ExtractBbox(void *pmodel, int sequence, float *mins, float *maxs)
{
	studiohdr_t *pstudiohdr = (studiohdr_t *)pmodel;

	if (!pstudiohdr)
		return 0;

	mstudioseqdesc_t *pseqdesc = (mstudioseqdesc_t *)((byte *)pstudiohdr + pstudiohdr->seqindex);

	mins[0] = pseqdesc[sequence].bbmin[0];
	mins[1] = pseqdesc[sequence].bbmin[1];
	mins[2] = pseqdesc[sequence].bbmin[2];
	maxs[0] = pseqdesc[sequence].bbmax[0];
	maxs[1] = pseqdesc[sequence].bbmax[1];
	maxs[2] = pseqdesc[sequence].bbmax[2];
	return 1;
}

int LookupActivity(void *pmodel, entvars_t *pev, int activity)
{
	studiohdr_t *pstudiohdr = (studiohdr_t *)pmodel;

	if (!pstudiohdr)
		return 0;

	mstudioseqdesc_t *pseqdesc = (mstudioseqdesc_t *)((byte *)pstudiohdr + pstudiohdr->seqindex);
	int weighttotal = 0;
	int seq = ACTIVITY_NOT_AVAILABLE;

	for (int i = 0; i < pstudiohdr->numseq; i++)
	{
		if (pseqdesc[i].activity == activity)
		{
			weighttotal += pseqdesc[i].actweight;

			if (!weighttotal || RANDOM_LONG(0, weighttotal - 1) < pseqdesc[i].actweight)
				seq = i;
		}
	}

	return seq;
}

int LookupActivityHeaviest(void *pmodel, entvars_t *pev, int activity)
{
	studiohdr_t *pstudiohdr = (studiohdr_t *)pmodel;

	if (!pstudiohdr)
		return 0;

	mstudioseqdesc_t *pseqdesc = (mstudioseqdesc_t *)((byte *)pstudiohdr + pstudiohdr->seqindex);
	int weight = 0;
	int seq = ACTIVITY_NOT_AVAILABLE;

	for (int i = 0; i < pstudiohdr->numseq; i++)
	{
		if (pseqdesc[i].activity == activity)
		{
			if (pseqdesc[i].actweight > weight)
			{
				weight = pseqdesc[i].actweight;
				seq = i;
			}
		}
	}

	return seq;
}

void GetEyePosition(void *pmodel, float *vecEyePosition)
{
	studiohdr_t *pstudiohdr = (studiohdr_t *)pmodel;

	if (!pstudiohdr)
	{
		ALERT(at_console, "GetEyePosition() Can't get pstudiohdr ptr!\n");
		return;
	}

	VectorCopy(pstudiohdr->eyeposition, vecEyePosition);
}

int LookupSequence(void *pmodel, const char *label)
{
	studiohdr_t *pstudiohdr = (studiohdr_t *)pmodel;

	if (!pstudiohdr)
		return 0;

	mstudioseqdesc_t *pseqdesc = (mstudioseqdesc_t *)((byte *)pstudiohdr + pstudiohdr->seqindex);

	for (int i = 0; i < pstudiohdr->numseq; i++)
	{
		if (!stricmp(pseqdesc[i].label, label))
			return i;
	}

	return -1;
}

int IsSoundEvent(int eventNumber)
{
	if (eventNumber == SCRIPT_EVENT_SOUND || eventNumber == SCRIPT_EVENT_SOUND_VOICE)
		return 1;

	return 0;
}

void SequencePrecache(void *pmodel, const char *pSequenceName)
{
	int index = LookupSequence(pmodel, pSequenceName);

	if (index >= 0)
	{
		studiohdr_t *pstudiohdr = (studiohdr_t *)pmodel;

		if (!pstudiohdr || index >= pstudiohdr->numseq)
			return;

		mstudioseqdesc_t *pseqdesc = (mstudioseqdesc_t *)((byte *)pstudiohdr + pstudiohdr->seqindex) + index;
		mstudioevent_t *pevent = (mstudioevent_t *)((byte *)pstudiohdr + pseqdesc->eventindex);

		for (int i = 0; i < pseqdesc->numevents; i++)
		{
			if (pevent[i].event >= EVENT_CLIENT)
				continue;

			if (IsSoundEvent(pevent[i].event))
			{
				if (!strlen(pevent[i].options))
					ALERT(at_error, "Bad sound event %d in sequence %s :: %s (sound is \"%s\")\n", pevent[i].event, pstudiohdr->name, pSequenceName, pevent[i].options);

				PRECACHE_SOUND((char *)(gpGlobals->pStringBase + ALLOC_STRING(pevent[i].options)));
			}
		}
	}
}

void GetSequenceInfo(void *pmodel, entvars_t *pev, float *pflFrameRate, float *pflGroundSpeed)
{
	studiohdr_t *pstudiohdr = (studiohdr_t *)pmodel;

	if (!pstudiohdr)
		return;

	if (pev->sequence >= pstudiohdr->numseq)
	{
		*pflFrameRate = 0;
		*pflGroundSpeed = 0;
		return;
	}

	mstudioseqdesc_t *pseqdesc = (mstudioseqdesc_t *)((byte *)pstudiohdr + pstudiohdr->seqindex) + (int)pev->sequence;

	if (pseqdesc->numframes > 1)
	{
		*pflFrameRate = 256 * pseqdesc->fps / (pseqdesc->numframes - 1);
		*pflGroundSpeed = sqrt(pseqdesc->linearmovement[0] * pseqdesc->linearmovement[0] + pseqdesc->linearmovement[1] * pseqdesc->linearmovement[1] + pseqdesc->linearmovement[2] * pseqdesc->linearmovement[2]);
		*pflGroundSpeed = *pflGroundSpeed * pseqdesc->fps / (pseqdesc->numframes - 1);
	}
	else
	{
		*pflFrameRate = 256;
		*pflGroundSpeed = 0;
	}
}

int GetSequenceFlags(void *pmodel, entvars_t *pev)
{
	studiohdr_t *pstudiohdr = (studiohdr_t *)pmodel;

	if (!pstudiohdr || pev->sequence >= pstudiohdr->numseq)
		return 0;

	mstudioseqdesc_t *pseqdesc = (mstudioseqdesc_t *)((byte *)pstudiohdr + pstudiohdr->seqindex) + (int)pev->sequence;
	return pseqdesc->flags;
}

int GetAnimationEvent(void *pmodel, entvars_t *pev, MonsterEvent_t *pMonsterEvent, float flStart, float flEnd, int index)
{
	studiohdr_t *pstudiohdr = (studiohdr_t *)pmodel;

	if (!pstudiohdr || pev->sequence >= pstudiohdr->numseq || !pMonsterEvent)
		return 0;

	mstudioseqdesc_t *pseqdesc = (mstudioseqdesc_t *)((byte *)pstudiohdr + pstudiohdr->seqindex) + (int)pev->sequence;
	mstudioevent_t *pevent = (mstudioevent_t *)((byte *)pstudiohdr + pseqdesc->eventindex);

	if (!pseqdesc->numevents || index > pseqdesc->numevents)
		return 0;

	if (pseqdesc->numframes > 1)
	{
		flStart *= (pseqdesc->numframes - 1) / 256;
		flEnd *= (pseqdesc->numframes - 1) / 256;
	}
	else
	{
		flStart = 0;
		flEnd = 1;
	}

	for (; index < pseqdesc->numevents; index++)
	{
		if (pevent[index].event >= EVENT_CLIENT)
			continue;

		if ((pevent[index].frame >= flStart && pevent[index].frame < flEnd) || ((pseqdesc->flags & STUDIO_LOOPING) && flEnd >= pseqdesc->numframes - 1 && pevent[index].frame < flEnd - pseqdesc->numframes + 1))
		{
			pMonsterEvent->event = pevent[index].event;
			pMonsterEvent->options = pevent[index].options;
			return index + 1;
		}
	}

	return 0;
}

float SetController(void *pmodel, entvars_t *pev, int iController, float flValue)
{
	studiohdr_t *pstudiohdr = (studiohdr_t *)pmodel;

	if (!pstudiohdr)
		return flValue;

	mstudiobonecontroller_t *pbonecontroller = (mstudiobonecontroller_t *)((byte *)pstudiohdr + pstudiohdr->bonecontrollerindex);
	int i = 0;

	for (; i < pstudiohdr->numbonecontrollers; i++, pbonecontroller++)
	{
		if (pbonecontroller->index == iController)
			break;
	}

	if (i >= pstudiohdr->numbonecontrollers)
		return flValue;

	if (pbonecontroller->type & (STUDIO_XR | STUDIO_YR | STUDIO_ZR))
	{
		if (pbonecontroller->end < pbonecontroller->start)
			flValue = -flValue;

		if (pbonecontroller->start + 359 >= pbonecontroller->end)
		{
			if (flValue > ((pbonecontroller->start + pbonecontroller->end) / 2) + 180)
				flValue = flValue - 360;

			if (flValue < ((pbonecontroller->start + pbonecontroller->end) / 2) - 180)
				flValue = flValue + 360;
		}
		else
		{
			if (flValue > 360)
				flValue = flValue - (int)(flValue / 360.0) * 360.0;
			else if (flValue < 0)
				flValue = flValue + (int)((flValue / -360.0) + 1) * 360.0;
		}
	}

	int setting = (int)(255 * (flValue - pbonecontroller->start) / (pbonecontroller->end - pbonecontroller->start));

	if (setting < 0)
		setting = 0;

	if (setting > 255)
		setting = 255;

	pev->controller[iController] = setting;
	return setting * (1.0 / 255) * (pbonecontroller->end - pbonecontroller->start) + pbonecontroller->start;
}

float SetBlending(void *pmodel, entvars_t *pev, int iBlender, float flValue)
{
	studiohdr_t *pstudiohdr = (studiohdr_t *)pmodel;

	if (!pstudiohdr)
		return flValue;

	mstudioseqdesc_t *pseqdesc = (mstudioseqdesc_t *)((byte *)pstudiohdr + pstudiohdr->seqindex) + (int)pev->sequence;

	if (pseqdesc->blendtype[iBlender] == 0)
		return flValue;

	if (pseqdesc->blendtype[iBlender] & (STUDIO_XR | STUDIO_YR | STUDIO_ZR))
	{
		if (pseqdesc->blendend[iBlender] < pseqdesc->blendstart[iBlender])
			flValue = -flValue;

		if (pseqdesc->blendstart[iBlender] + 359 >= pseqdesc->blendend[iBlender])
		{
			if (flValue > ((pseqdesc->blendstart[iBlender] + pseqdesc->blendend[iBlender]) / 2) + 180)
				flValue = flValue - 360;

			if (flValue < ((pseqdesc->blendstart[iBlender] + pseqdesc->blendend[iBlender]) / 2) - 180)
				flValue = flValue + 360;
		}
	}

	int setting = (int)(255 * (flValue - pseqdesc->blendstart[iBlender]) / (pseqdesc->blendend[iBlender] - pseqdesc->blendstart[iBlender]));

	if (setting < 0)
		setting = 0;

	if (setting > 255)
		setting = 255;

	pev->blending[iBlender] = setting;
	return setting * (1.0 / 255) * (pseqdesc->blendend[iBlender] - pseqdesc->blendstart[iBlender]) + pseqdesc->blendstart[iBlender];
}

int FindTransition(void *pmodel, int iEndingAnim, int iGoalAnim, int *piDir)
{
	studiohdr_t *pstudiohdr = (studiohdr_t *)pmodel;

	if (!pstudiohdr)
		return iGoalAnim;

	mstudioseqdesc_t *pseqdesc = (mstudioseqdesc_t *)((byte *)pstudiohdr + pstudiohdr->seqindex);

	if (!pseqdesc[iEndingAnim].entrynode || pseqdesc[iGoalAnim].entrynode)
		return iGoalAnim;

	int iEndNode = (*piDir > 0) ? pseqdesc[iEndingAnim].exitnode : pseqdesc[iEndingAnim].entrynode;

	if (iEndNode == pseqdesc[iGoalAnim].entrynode)
	{
		*piDir = 1;
		return iGoalAnim;
	}

	byte *pTransition = ((byte *)pstudiohdr + pstudiohdr->transitionindex);
	int iInternNode = pTransition[(iEndNode - 1) * pstudiohdr->numtransitions + (pseqdesc[iGoalAnim].entrynode - 1)];

	if (iInternNode == 0)
		return iGoalAnim;

	for (int i = 0; i < pstudiohdr->numseq; i++)
	{
		if (pseqdesc[i].entrynode == iEndNode && pseqdesc[i].exitnode == iInternNode)
		{
			*piDir = 1;
			return i;
		}

		if (pseqdesc[i].nodeflags)
		{
			if (pseqdesc[i].exitnode == iEndNode && pseqdesc[i].entrynode == iInternNode)
			{
				*piDir = -1;
				return i;
			}
		}
	}

	ALERT(at_console, "error in transition graph");
	return iGoalAnim;
}

void SetBodygroup(void *pmodel, entvars_t *pev, int iGroup, int iValue)
{
	studiohdr_t *pstudiohdr = (studiohdr_t *)pmodel;

	if (!pstudiohdr)
		return;

	if (iGroup > pstudiohdr->numbodyparts)
		return;

	mstudiobodyparts_t *pbodypart = (mstudiobodyparts_t *)((byte *)pstudiohdr + pstudiohdr->bodypartindex) + iGroup;

	if (iValue >= pbodypart->nummodels)
		return;

	int iCurrent = (pev->body / pbodypart->base) % pbodypart->nummodels;
	pev->body = (pev->body - (iCurrent * pbodypart->base) + (iValue * pbodypart->base));
}

int GetBodygroup(void *pmodel, entvars_t *pev, int iGroup)
{
	studiohdr_t *pstudiohdr = (studiohdr_t *)pmodel;

	if (!pstudiohdr)
		return 0;

	if (iGroup > pstudiohdr->numbodyparts)
		return 0;

	mstudiobodyparts_t *pbodypart = (mstudiobodyparts_t *)((byte *)pstudiohdr + pstudiohdr->bodypartindex) + iGroup;

	if (pbodypart->nummodels <= 1)
		return 0;

	int iCurrent = (pev->body / pbodypart->base) % pbodypart->nummodels;
	return iCurrent;
}

studiohdr_t *g_pstudiohdr;
void SV_StudioSetupBones(struct model_s *pModel, float frame, int sequence, const vec3_t angles, const vec3_t origin, const byte *pcontroller, const byte *pblending, int iBone, const edict_t *pEdict);

sv_blending_interface_t svBlending =
{
	SV_BLENDING_INTERFACE_VERSION,
	SV_StudioSetupBones
};

server_studio_api_t IEngineStudio;

float (*g_pRotationMatrix)[3][4];
float (*g_pBoneTransform)[MAXSTUDIOBONES][3][4];

extern "C" EXPORT int Server_GetBlendingInterface(int version, sv_blending_interface_t **pinterface, server_studio_api_t *pstudio, float ***rotationmatrix, float ****bonetransform)
{
	if (version != SV_BLENDING_INTERFACE_VERSION)
		return 0;

	*pinterface = &svBlending;
	IEngineStudio.Mem_Calloc = pstudio->Mem_Calloc;
	IEngineStudio.Cache_Check = pstudio->Cache_Check;
	IEngineStudio.LoadCacheFile = pstudio->LoadCacheFile;
	IEngineStudio.Mod_Extradata = pstudio->Mod_Extradata;
	g_pRotationMatrix = (float (*)[3][4])rotationmatrix;
	g_pBoneTransform = (float (*)[MAXSTUDIOBONES][3][4])bonetransform;
	return 1;
}

void AngleQuaternion(const vec3_t angles, vec4_t quaternion)
{
	float angle;
	float sr, sp, sy, cr, cp, cy;

	angle = angles[2] * 0.5;
	sy = sin(angle);
	cy = cos(angle);
	angle = angles[1] * 0.5;
	sp = sin(angle);
	cp = cos(angle);
	angle = angles[0] * 0.5;
	sr = sin(angle);
	cr = cos(angle);

	quaternion[0] = sr * cp * cy - cr * sp * sy;
	quaternion[1] = cr * sp * cy + sr * cp * sy;
	quaternion[2] = cr * cp * sy - sr * sp * cy;
	quaternion[3] = cr * cp * cy + sr * sp * sy;
}

void QuaternionSlerp(const vec4_t p, vec4_t q, float t, vec4_t qt)
{
	int i;
	float omega, cosom, sinom, sclp, sclq;
	float a = 0;
	float b = 0;

	for (i = 0; i < 4; i++)
	{
		a += (p[i]-q[i]) * (p[i]-q[i]);
		b += (p[i]+q[i]) * (p[i]+q[i]);
	}

	if (a > b)
	{
		for (i = 0; i < 4; i++)
			q[i] = -q[i];
	}

	cosom = p[0] * q[0] + p[1] * q[1] + p[2] * q[2] + p[3] * q[3];

	if ((1.0 + cosom) > 0.00000001)
	{
		if ((1.0 - cosom) > 0.00000001)
		{
			omega = acos(cosom);
			sinom = sin(omega);
			sclp = sin((1.0 - t) * omega) / sinom;
			sclq = sin(t * omega) / sinom;
		}
		else
		{
			sclp = 1.0 - t;
			sclq = t;
		}

		for (i = 0; i < 4; i++)
			qt[i] = sclp * p[i] + sclq * q[i];
	}
	else
	{
		qt[0] = -p[1];
		qt[1] = p[0];
		qt[2] = -p[3];
		qt[3] = p[2];
		sclp = sin((1.0 - t) * 0.5 * M_PI);
		sclq = sin(t * 0.5 * M_PI);

		for (i = 0; i < 3; i++)
			qt[i] = sclp * p[i] + sclq * qt[i];
	}
}

void QuaternionMatrix(const vec4_t quaternion, float (*matrix)[4])
{
	matrix[0][0] = 1.0 - 2.0 * quaternion[1] * quaternion[1] - 2.0 * quaternion[2] * quaternion[2];
	matrix[1][0] = 2.0 * quaternion[0] * quaternion[1] + 2.0 * quaternion[3] * quaternion[2];
	matrix[2][0] = 2.0 * quaternion[0] * quaternion[2] - 2.0 * quaternion[3] * quaternion[1];
	matrix[0][1] = 2.0 * quaternion[0] * quaternion[1] - 2.0 * quaternion[3] * quaternion[2];
	matrix[1][1] = 1.0 - 2.0 * quaternion[0] * quaternion[0] - 2.0 * quaternion[2] * quaternion[2];
	matrix[2][1] = 2.0 * quaternion[1] * quaternion[2] + 2.0 * quaternion[3] * quaternion[0];
	matrix[0][2] = 2.0 * quaternion[0] * quaternion[2] + 2.0 * quaternion[3] * quaternion[1];
	matrix[1][2] = 2.0 * quaternion[1] * quaternion[2] - 2.0 * quaternion[3] * quaternion[0];
	matrix[2][2] = 1.0 - 2.0 * quaternion[0] * quaternion[0] - 2.0 * quaternion[1] * quaternion[1];
}

mstudioanim_t *StudioGetAnim(model_t *model, mstudioseqdesc_t *pseqdesc)
{
	mstudioseqgroup_t *pseqgroup = (mstudioseqgroup_t *)((byte *)g_pstudiohdr + g_pstudiohdr->seqgroupindex) + pseqdesc->seqgroup;

	if (!pseqdesc->seqgroup)
		return (mstudioanim_t *)((byte *)g_pstudiohdr + pseqgroup->data + pseqdesc->animindex);

	cache_user_t *paSequences = (cache_user_t *)model->submodels;

	if (!paSequences)
	{
		paSequences = (cache_user_t *)IEngineStudio.Mem_Calloc(16, sizeof(cache_user_t));
		model->submodels = (dmodel_t *)paSequences;
	}

	if (!IEngineStudio.Cache_Check(&paSequences[pseqdesc->seqgroup]))
		IEngineStudio.LoadCacheFile(pseqgroup->name, &paSequences[pseqdesc->seqgroup]);

	return (mstudioanim_t *)((byte *)paSequences[pseqdesc->seqgroup].data + pseqdesc->animindex);
}

void StudioCalcBoneAdj(float dadt, float *adj, const byte *pcontroller1, const byte *pcontroller2, byte mouthopen)
{
	mstudiobonecontroller_t *pbonecontroller = (mstudiobonecontroller_t *)((byte *)g_pstudiohdr + g_pstudiohdr->bonecontrollerindex);

	for (int i = 0; i < g_pstudiohdr->numbonecontrollers; i++)
	{
		float value;
		int index = pbonecontroller[i].index;

		if (index <= 3)
		{
			if (pbonecontroller[i].type & STUDIO_RLOOP)
			{
				if (abs(pcontroller1[index] - pcontroller2[index]) > 128)
					value = ((((pcontroller1[i] + 128) % 256) * dadt) + (((pcontroller2[i] + 128) % 256) * (1 - dadt)) - 128) * (360.0 / 256) + pbonecontroller[i].start;
				else
					value = ((1.0 - dadt) * pcontroller1[index] + pcontroller2[index] * dadt) * (360.0 / 256) + pbonecontroller[i].start;
			}
			else
			{
				value = (pcontroller1[index] * dadt + pcontroller2[index] * (1.0 - dadt)) / 255.0;

				if (value < 0)
					value = 0;

				if (value > 1)
					value = 1;

				value = (1.0 - value) * pbonecontroller[i].start + value * pbonecontroller[i].end;
			}
		}
		else
		{
			value = mouthopen / 64;

			if (value > 1)
				value = 1;

			value = (1.0 - value) * pbonecontroller[i].start + value * pbonecontroller[i].end;
		}

		switch (pbonecontroller[i].type & STUDIO_TYPES)
		{
			case STUDIO_XR:
			case STUDIO_YR:
			case STUDIO_ZR: adj[i] = value * (M_PI / 180.0); break;
			case STUDIO_X:
			case STUDIO_Y:
			case STUDIO_Z: adj[i] = value; break;
		}
	}
}

void StudioCalcBoneQuaterion(int frame, float s, mstudiobone_t *pbone, mstudioanim_t *panim, float *adj, float *q)
{
	vec4_t q1, q2;
	vec3_t angle1, angle2;

	for (int i = 0; i < 3; i++)
	{
		if (panim->offset[i + 3])
		{
			mstudioanimvalue_t *panimvalue = (mstudioanimvalue_t *)((byte *)panim + panim->offset[i + 3]);
			int j = (panimvalue->num.total < panimvalue->num.valid) ? 0 : frame;

			while (panimvalue->num.total <= j)
			{
				j -= panimvalue->num.total;
				panimvalue += panimvalue->num.valid + 1;

				if (panimvalue->num.total < panimvalue->num.valid)
					j = 0;
			}

			if (panimvalue->num.valid > j)
			{
				angle1[i] = panimvalue[j + 1].value;

				if (panimvalue->num.valid > j + 1)
					angle2[i] = panimvalue[j + 2].value;
				else if (panimvalue->num.total > j + 1)
					angle2[i] = angle1[i];
				else
					angle2[i] = panimvalue[panimvalue->num.valid + 2].value;
			}
			else
			{
				angle1[i] = panimvalue[panimvalue->num.valid].value;

				if (panimvalue->num.total > j + 1)
					angle2[i] = angle1[i];
				else
					angle2[i] = panimvalue[panimvalue->num.valid + 2].value;
			}

			angle1[i] = pbone->value[i + 3] + angle1[i] * pbone->scale[i + 3];
			angle2[i] = pbone->value[i + 3] + angle2[i] * pbone->scale[i + 3];
		}
		else
			angle2[i] = angle1[i] = pbone->value[i + 3];

		if (pbone->bonecontroller[i + 3] != -1)
		{
			angle1[i] += adj[pbone->bonecontroller[i + 3]];
			angle2[i] += adj[pbone->bonecontroller[i + 3]];
		}
	}

	if (!VectorCompare(angle1, angle2))
	{
		AngleQuaternion(angle1, q1);
		AngleQuaternion(angle2, q2);
		QuaternionSlerp(q1, q2, s, q);
	}
	else
		AngleQuaternion(angle1, q);
}

void StudioCalcBonePosition(int frame, float s, mstudiobone_t *pbone, mstudioanim_t *panim, float *adj, float *pos)
{
	for (int i = 0; i < 3; i++)
	{
		pos[i] = pbone->value[i];

		if (panim->offset[i] != 0)
		{
			mstudioanimvalue_t *panimvalue = (mstudioanimvalue_t *)((byte *)panim + panim->offset[i]);
			int j = (panimvalue->num.total < panimvalue->num.valid) ? 0 : frame;

			while (panimvalue->num.total <= j)
			{
				j -= panimvalue->num.total;
				panimvalue += panimvalue->num.valid + 1;

				if (panimvalue->num.total < panimvalue->num.valid)
					j = 0;
			}

			if (panimvalue->num.valid > j)
			{
				if (panimvalue->num.valid > j + 1)
					pos[i] += (panimvalue[j + 1].value * (1 - s) + s * panimvalue[j + 2].value) * pbone->scale[i];
				else
					pos[i] += panimvalue[j + 1].value * pbone->scale[i];
			}
			else
			{
				if (panimvalue->num.total <= j + 1)
					pos[i] += (panimvalue[panimvalue->num.valid].value * (1 - s) + s * panimvalue[panimvalue->num.valid + 2].value) * pbone->scale[i];
				else
					pos[i] += panimvalue[panimvalue->num.valid].value * pbone->scale[i];
			}
		}

		if (pbone->bonecontroller[i] != -1 && adj)
			pos[i] += adj[pbone->bonecontroller[i]];
	}
}

void StudioSlerpBones(vec4_t q1[], float pos1[][3], vec4_t q2[], float pos2[][3], float s)
{
	if (s < 0)
		s = 0;
	else if (s > 1)
		s = 1;

	float s1 = 1 - s;

	for (int i = 0; i < g_pstudiohdr->numbones; i++)
	{
		vec4_t q3;
		QuaternionSlerp(q1[i], q2[i], s, q3);

		q1[i][0] = q3[0];
		q1[i][1] = q3[1];
		q1[i][2] = q3[2];
		q1[i][3] = q3[3];
		pos1[i][0] = pos1[i][0] * s1 + pos2[i][0] * s;
		pos1[i][1] = pos1[i][1] * s1 + pos2[i][1] * s;
		pos1[i][2] = pos1[i][2] * s1 + pos2[i][2] * s;
	}
}

float StudioEstimateFrame(float frame, mstudioseqdesc_t *pseqdesc)
{
	if (pseqdesc->numframes <= 1)
		return 0;

	return frame * (pseqdesc->numframes - 1) / 256;
}

mstudioanim_t *LookupAnimation(model_t *model, mstudioseqdesc_t *pseqdesc, int index)
{
	mstudioanim_t *panim = NULL;

	panim = StudioGetAnim(model, pseqdesc);

	assert(panim);

	if (index < 0)
		return panim;

	if (index > (pseqdesc->numblends - 1))
		return panim;

	panim += index * g_pstudiohdr->numbones;
	return panim;
}

void ConcatTransforms(const float in1[3][4], const float in2[3][4], float out[3][4])
{
	out[0][0] = in1[0][0] * in2[0][0] + in1[0][1] * in2[1][0] + in1[0][2] * in2[2][0];
	out[0][1] = in1[0][0] * in2[0][1] + in1[0][1] * in2[1][1] + in1[0][2] * in2[2][1];
	out[0][2] = in1[0][0] * in2[0][2] + in1[0][1] * in2[1][2] + in1[0][2] * in2[2][2];
	out[0][3] = in1[0][0] * in2[0][3] + in1[0][1] * in2[1][3] + in1[0][2] * in2[2][3] + in1[0][3];
	out[1][0] = in1[1][0] * in2[0][0] + in1[1][1] * in2[1][0] + in1[1][2] * in2[2][0];
	out[1][1] = in1[1][0] * in2[0][1] + in1[1][1] * in2[1][1] + in1[1][2] * in2[2][1];
	out[1][2] = in1[1][0] * in2[0][2] + in1[1][1] * in2[1][2] + in1[1][2] * in2[2][2];
	out[1][3] = in1[1][0] * in2[0][3] + in1[1][1] * in2[1][3] + in1[1][2] * in2[2][3] + in1[1][3];
	out[2][0] = in1[2][0] * in2[0][0] + in1[2][1] * in2[1][0] + in1[2][2] * in2[2][0];
	out[2][1] = in1[2][0] * in2[0][1] + in1[2][1] * in2[1][1] + in1[2][2] * in2[2][1];
	out[2][2] = in1[2][0] * in2[0][2] + in1[2][1] * in2[1][2] + in1[2][2] * in2[2][2];
	out[2][3] = in1[2][0] * in2[0][3] + in1[2][1] * in2[1][3] + in1[2][2] * in2[2][3] + in1[2][3];
}

float GetPlayerPitch(const edict_t *pent);
float GetPlayerYaw(const edict_t *pent);
int GetPlayerGaitsequence(const edict_t *pent);
float UTIL_GetPlayerGaitYaw(int playerIndex);

void SV_StudioSetupBones(struct model_s *pModel, float frame, int sequence, const vec3_t angles, const vec3_t origin, const byte *pcontroller, const byte *pblending, int iBone, const edict_t *pEdict)
{
	int i, j;
	float f;
	int subframe;
	float adj[MAXSTUDIOCONTROLLERS];
	mstudiobone_t *pbones;
	mstudioseqdesc_t *pseqdesc;
	mstudioanim_t *panim;
	float bonematrix[3][4];
	int chain[MAXSTUDIOBONES];
	int chainlength;
	vec3_t temp_angles;

	static float pos[MAXSTUDIOBONES][3], pos2[MAXSTUDIOBONES][3], pos3[MAXSTUDIOBONES][3], pos4[MAXSTUDIOBONES][3];
	static vec4_t q[MAXSTUDIOBONES], q2[MAXSTUDIOBONES], q3[MAXSTUDIOBONES], q4[MAXSTUDIOBONES];

	g_pstudiohdr = (studiohdr_t *)IEngineStudio.Mod_Extradata(pModel);

	if (sequence < 0 || sequence >= g_pstudiohdr->numseq)
		sequence = 0;

	pbones = (mstudiobone_t *)((byte *)g_pstudiohdr + g_pstudiohdr->boneindex);
	pseqdesc = (mstudioseqdesc_t *)((byte *)g_pstudiohdr + g_pstudiohdr->seqindex) + sequence;
	panim = StudioGetAnim(pModel, pseqdesc);

	if (iBone < -1 || iBone >= g_pstudiohdr->numbones)
		iBone = 0;

	if (iBone == -1)
	{
		chainlength = g_pstudiohdr->numbones;

		for (i = 0; i < chainlength; i++)
			chain[(chainlength - i) - 1] = i;
	}
	else
	{
		chainlength = 0;

		for (i = iBone; i != -1; i = pbones[i].parent)
			chain[chainlength++] = i;
	}

	f = StudioEstimateFrame(frame, pseqdesc);
	subframe = f;
	f -= subframe;
	StudioCalcBoneAdj(0, adj, pcontroller, pcontroller, 0);

	for (i = chainlength - 1; i >= 0; i--)
	{
		j = chain[i];
		StudioCalcBoneQuaterion(subframe, f, &pbones[j], &panim[j], adj, q[j]);
		StudioCalcBonePosition(subframe, f, &pbones[j], &panim[j], adj, pos[j]);
	}

	if (pseqdesc->numblends == 9)
	{
		float s, t;

		s = GetPlayerYaw(pEdict);
		t = GetPlayerPitch(pEdict);

		if (s <= 127.0)
		{
			s = (s * 2.0);

			if (t <= 127.0)
			{
				t = (t * 2.0);

				for (i = chainlength - 1; i >= 0; i--)
				{
					j = chain[i];
					StudioCalcBoneQuaterion(subframe, f, &pbones[j], &panim[j], adj, q[j]);
					StudioCalcBonePosition(subframe, f, &pbones[j], &panim[j], adj, pos[j]);
				}

				panim = LookupAnimation(pModel, pseqdesc, 1);

				for (i = chainlength - 1; i >= 0; i--)
				{
					j = chain[i];
					StudioCalcBoneQuaterion(subframe, f, &pbones[j], &panim[j], adj, q2[j]);
					StudioCalcBonePosition(subframe, f, &pbones[j], &panim[j], adj, pos2[j]);
				}

				panim = LookupAnimation(pModel, pseqdesc, 3);

				for (i = chainlength - 1; i >= 0; i--)
				{
					j = chain[i];
					StudioCalcBoneQuaterion(subframe, f, &pbones[j], &panim[j], adj, q3[j]);
					StudioCalcBonePosition(subframe, f, &pbones[j], &panim[j], adj, pos3[j]);
				}

				panim = LookupAnimation(pModel, pseqdesc, 4);

				for (i = chainlength - 1; i >= 0; i--)
				{
					j = chain[i];
					StudioCalcBoneQuaterion(subframe, f, &pbones[j], &panim[j], adj, q4[j]);
					StudioCalcBonePosition(subframe, f, &pbones[j], &panim[j], adj, pos4[j]);
				}
			}
			else
			{
				t = 2.0 * (t - 127.0);

				panim = LookupAnimation(pModel, pseqdesc, 3);

				for (i = chainlength - 1; i >= 0; i--)
				{
					j = chain[i];
					StudioCalcBoneQuaterion(subframe, f, &pbones[j], &panim[j], adj, q[j]);
					StudioCalcBonePosition(subframe, f, &pbones[j], &panim[j], adj, pos[j]);
				}

				panim = LookupAnimation(pModel, pseqdesc, 4);

				for (i = chainlength - 1; i >= 0; i--)
				{
					j = chain[i];
					StudioCalcBoneQuaterion(subframe, f, &pbones[j], &panim[j], adj, q2[j]);
					StudioCalcBonePosition(subframe, f, &pbones[j], &panim[j], adj, pos2[j]);
				}

				panim = LookupAnimation(pModel, pseqdesc, 6);

				for (i = chainlength - 1; i >= 0; i--)
				{
					j = chain[i];
					StudioCalcBoneQuaterion(subframe, f, &pbones[j], &panim[j], adj, q3[j]);
					StudioCalcBonePosition(subframe, f, &pbones[j], &panim[j], adj, pos3[j]);
				}

				panim = LookupAnimation(pModel, pseqdesc, 7);

				for (i = chainlength - 1; i >= 0; i--)
				{
					j = chain[i];
					StudioCalcBoneQuaterion(subframe, f, &pbones[j], &panim[j], adj, q4[j]);
					StudioCalcBonePosition(subframe, f, &pbones[j], &panim[j], adj, pos4[j]);
				}
			}
		}
		else
		{
			s = 2.0 * (s - 127.0);

			if (t <= 127.0)
			{
				t = (t * 2.0);

				panim = LookupAnimation(pModel, pseqdesc, 1);

				for (i = chainlength - 1; i >= 0; i--)
				{
					j = chain[i];
					StudioCalcBoneQuaterion(subframe, f, &pbones[j], &panim[j], adj, q[j]);
					StudioCalcBonePosition(subframe, f, &pbones[j], &panim[j], adj, pos[j]);
				}

				panim = LookupAnimation(pModel, pseqdesc, 2);

				for (i = chainlength - 1; i >= 0; i--)
				{
					j = chain[i];
					StudioCalcBoneQuaterion(subframe, f, &pbones[j], &panim[j], adj, q2[j]);
					StudioCalcBonePosition(subframe, f, &pbones[j], &panim[j], adj, pos2[j]);
				}

				panim = LookupAnimation(pModel, pseqdesc, 4);

				for (i = chainlength - 1; i >= 0; i--)
				{
					j = chain[i];
					StudioCalcBoneQuaterion(subframe, f, &pbones[j], &panim[j], adj, q3[j]);
					StudioCalcBonePosition(subframe, f, &pbones[j], &panim[j], adj, pos3[j]);
				}

				panim = LookupAnimation(pModel, pseqdesc, 5);

				for (i = chainlength - 1; i >= 0; i--)
				{
					j = chain[i];
					StudioCalcBoneQuaterion(subframe, f, &pbones[j], &panim[j], adj, q4[j]);
					StudioCalcBonePosition(subframe, f, &pbones[j], &panim[j], adj, pos4[j]);
				}
			}
			else
			{
				t = 2.0 * (t - 127.0);

				panim = LookupAnimation(pModel, pseqdesc, 4);

				for (i = chainlength - 1; i >= 0; i--)
				{
					j = chain[i];
					StudioCalcBoneQuaterion(subframe, f, &pbones[j], &panim[j], adj, q[j]);
					StudioCalcBonePosition(subframe, f, &pbones[j], &panim[j], adj, pos[j]);
				}

				panim = LookupAnimation(pModel, pseqdesc, 5);

				for (i = chainlength - 1; i >= 0; i--)
				{
					j = chain[i];
					StudioCalcBoneQuaterion(subframe, f, &pbones[j], &panim[j], adj, q2[j]);
					StudioCalcBonePosition(subframe, f, &pbones[j], &panim[j], adj, pos2[j]);
				}

				panim = LookupAnimation(pModel, pseqdesc, 7);

				for (i = chainlength - 1; i >= 0; i--)
				{
					j = chain[i];
					StudioCalcBoneQuaterion(subframe, f, &pbones[j], &panim[j], adj, q3[j]);
					StudioCalcBonePosition(subframe, f, &pbones[j], &panim[j], adj, pos3[j]);
				}

				panim = LookupAnimation(pModel, pseqdesc, 8);

				for (i = chainlength - 1; i >= 0; i--)
				{
					j = chain[i];
					StudioCalcBoneQuaterion(subframe, f, &pbones[j], &panim[j], adj, q4[j]);
					StudioCalcBonePosition(subframe, f, &pbones[j], &panim[j], adj, pos4[j]);
				}
			}
		}

		s /= 255.0;
		t /= 255.0;

		StudioSlerpBones(q, pos, q2, pos2, s);
		StudioSlerpBones(q3, pos3, q4, pos4, s);
		StudioSlerpBones(q, pos, q3, pos3, t);
	}
	else if (pseqdesc->numblends > 1)
	{
		float s;

		pseqdesc = (mstudioseqdesc_t *)((byte *)g_pstudiohdr + g_pstudiohdr->seqindex) + sequence;
		panim = StudioGetAnim(pModel, pseqdesc);
		panim += g_pstudiohdr->numbones * 1;

		for (i = chainlength - 1; i >= 0; i--)
		{
			j = chain[i];
			StudioCalcBoneQuaterion(subframe, f, &pbones[j], &panim[j], adj, q2[j]);
			StudioCalcBonePosition(subframe, f, &pbones[j], &panim[j], adj, pos2[j]);
		}

		s = (float)pblending[0] / 255.0;
		StudioSlerpBones(q, pos, q2, pos2, s);
	}

	if (pseqdesc->numblends == 9 && sequence < ANIM_FIRST_DEATH_SEQUENCE && sequence != ANIM_SWIM_1 && sequence != ANIM_SWIM_2)
	{
		int copy = 0;
		int gaitsequence = GetPlayerGaitsequence(pEdict);

		if (gaitsequence >= g_pstudiohdr->numseq)
			gaitsequence = 0;

		if (gaitsequence < 0)
			gaitsequence = 0;

		pseqdesc = (mstudioseqdesc_t *)((byte *)g_pstudiohdr + g_pstudiohdr->seqindex) + gaitsequence;
		panim = StudioGetAnim(pModel, pseqdesc);

		for (i = chainlength - 1; i >= 0; i--)
		{
			j = chain[i];
			StudioCalcBoneQuaterion(0, 0, &pbones[j], &panim[j], adj, q2[j]);
			StudioCalcBonePosition(0, 0, &pbones[j], &panim[j], adj, pos2[j]);
		}

		for (i = 0; i < g_pstudiohdr->numbones; i++)
		{
			if (strcmp(pbones[i].name, "Bip01 Spine") == 0)
				copy = 0;
			else if (strcmp(pbones[pbones[i].parent].name, "Bip01 Pelvis") == 0)
				copy = 1;

			if (copy)
			{
				memcpy(pos[i], pos2[i], sizeof(pos[i]));
				memcpy(q[i], q2[i], sizeof(q[i]));
			}
		}
	}

	VectorCopy(angles, temp_angles);

	if (pEdict)
	{
		temp_angles[1] = UTIL_GetPlayerGaitYaw(g_engfuncs.pfnIndexOfEdict(pEdict));

		if (temp_angles[1] < 0)
			temp_angles[1] += 360;
	}

	AngleMatrix(temp_angles, (*g_pRotationMatrix));

	(*g_pRotationMatrix)[0][3] = origin[0];
	(*g_pRotationMatrix)[1][3] = origin[1];
	(*g_pRotationMatrix)[2][3] = origin[2];

	for (i = chainlength - 1; i >= 0; i--)
	{
		j = chain[i];
		QuaternionMatrix(q[j], bonematrix);

		bonematrix[0][3] = pos[j][0];
		bonematrix[1][3] = pos[j][1];
		bonematrix[2][3] = pos[j][2];

		if (pbones[j].parent == -1)
			ConcatTransforms((*g_pRotationMatrix), bonematrix, (*g_pBoneTransform)[j]);
		else
			ConcatTransforms((*g_pBoneTransform)[pbones[j].parent], bonematrix, (*g_pBoneTransform)[j]);
	}
}