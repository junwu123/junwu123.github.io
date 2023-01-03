/*************************************************************************************

CUDAFEM

Author: Christian Dick

Copyright (c) Christian Dick

mailto:chr.dick@googlemail.com

This source code is property of Christian Dick. All rights reserved.
Unauthorized use prohibited.

*************************************************************************************/


//#define SOFT_SHADOWS


SamplerState samLinearBorder
{
	Filter = MIN_MAG_LINEAR_MIP_POINT;
	AddressU = BORDER;
	AddressV = BORDER;
	BorderColor = float4(1.0f, 0.0f, 0.0f, 0.0f);
};

RasterizerState rsCullNone
{
	FrontCounterClockwise = TRUE;
	CullMode = NONE;
};

RasterizerState rsCullFront
{
	FrontCounterClockwise = TRUE;
	CullMode = FRONT;
};

RasterizerState rsCullBack
{
	FrontCounterClockwise = TRUE;
	CullMode = BACK;
};

RasterizerState rsCullNoneWireFrame
{
	FillMode = WIREFRAME;
	FrontCounterClockwise = TRUE;
	CullMode = NONE;
};

RasterizerState rsCullNonePolygonOffset
{
	FrontCounterClockwise = TRUE;
	CullMode = NONE;
    DepthBias = 1.1f;
    SlopeScaledDepthBias = 4.0f;
};

RasterizerState rsDefault
{
};

DepthStencilState dsEnableDepth
{
	DepthEnable = TRUE;
};

DepthStencilState dsDisableDepth
{
	DepthEnable = FALSE;
};

BlendState bsDefault
{
};

BlendState bsAlpha
{
	BlendEnable[0] = TRUE;
	SrcBlend = ONE;
	DestBlend = INV_SRC_ALPHA;
	BlendOp = ADD;
};


cbuffer cbRarely
{
	float4x4	g_mProjection;
	float3		g_vLightPosCS;
	float		g_fSphereRadius; // Radius of spheres used to render surface vertices
};


cbuffer cbOncePerObject
{
	float4x4	g_mWorldView;
	float4x4	g_mTransInvWorldView;
	float3		g_vAmbientColor;
	float3		g_vDiffuseColor;
	float3		g_vSpecularColor;
};


cbuffer cbOncePerObjectWithFog
{
	float4x4	g_mWorld;		// Only used in vsRenderMeshWithFog and vsRenderMeshWithShadowsAndClipping
	float3		g_vSceneCenter;
	float3		g_vFogColor;
	float		g_fFogExponent;
	float		g_fFogDensity;
};


cbuffer cbOncePerObjectIntoShadowMap
{
	float4x4	g_mWorldLightSMProjection;
}


cbuffer cbOncePerObjectWithCutSurfaces
{
	float3		g_vAmbientColorCutSurfaces;
	float3		g_vDiffuseColorCutSurfaces;
	float3		g_vSpecularColorCutSurfaces;
}


cbuffer cbOncePerObjectWithShadow
{
	float4x4	g_mShadowMap;	// The consistent name would be "mInvViewLightProjectionShadowMap"
	float4x4	g_mInvViewLight;
}


cbuffer cbOncePerObjectWithClipping
{
	float		g_fClippingY;
}


Texture2D<float>	g_txShadowMap;



float4 Phong(float3 vPosCS, float3 vNormalCS, float3 vLightPosCS, float3 vAmbientColor, float3 vDiffuseColor, float3 vSpecularColor, float fSpecularExponent)
{
	float3 vLight = normalize(vLightPosCS - vPosCS);
	float3 vView = normalize(-vPosCS);
	float vl = dot(vNormalCS, vLight);
	float3 vReflect = (2.0f * vl) * vNormalCS - vLight;
	
	return float4(
		vAmbientColor
		+ vDiffuseColor * saturate(vl)
		+ vSpecularColor * pow(saturate(dot(vReflect, vNormalCS)), fSpecularExponent)
		, 1.0f);
}


void vsRenderMeshWithFaceNormals(
	float3 vPos : POSITION,
	out float4 vPosOut : SV_POSITION,
	out float3 vPosCSOut : POS_CS
	)
{
	vPosOut = mul(g_mWorldView, float4(vPos, 1.0f));
	vPosCSOut = vPosOut.xyz;
	vPosOut = mul(g_mProjection, vPosOut);
}


struct GSRenderMeshWithFaceNormalsInput
{
	float4 vPos : SV_POSITION;
	float3 vPosCS : POS_CS;
};

struct GSRenderMeshWithFaceNormalsOutput
{
	float4 vPos : SV_POSITION;
	float3 vPosCS : POS_CS;
	float3 vNormalCS : NORMAL_CS;
};


[maxvertexcount(3)]
void gsRenderMeshWithFaceNormals(
	triangle GSRenderMeshWithFaceNormalsInput input[3],
	inout TriangleStream<GSRenderMeshWithFaceNormalsOutput> stream
	)
{
	float3 vNormalCS = normalize( cross(input[1].vPosCS - input[0].vPosCS, input[2].vPosCS - input[0].vPosCS) );
	GSRenderMeshWithFaceNormalsOutput output;
	output.vPos = input[0].vPos;
	output.vPosCS = input[0].vPosCS;
	output.vNormalCS = vNormalCS;
	stream.Append(output);
	output.vPos = input[1].vPos;
	output.vPosCS = input[1].vPosCS;
	stream.Append(output);
	output.vPos = input[2].vPos;
	output.vPosCS = input[2].vPosCS;
	stream.Append(output);
}


void vsRenderMesh(
	float3 vPos : POSITION,
	float3 vNormal : NORMAL,
	out float4 vPosOut : SV_POSITION,
	out float3 vPosCSOut : POS_CS,
	out float3 vNormalCSOut : NORMAL_CS
	)
{
	vPosOut = mul(g_mWorldView, float4(vPos, 1.0f));
	vPosCSOut = vPosOut.xyz;
	vPosOut = mul(g_mProjection, vPosOut);
	vNormalCSOut = normalize( mul(g_mTransInvWorldView, float4(vNormal, 0.0f)).xyz );
}


void psRenderMesh(
	float4 vPos : SV_POSITION,
	float3 vPosCS : POS_CS,
	float3 vNormalCS : NORMAL_CS,
	out float4 vColorOut : SV_TARGET,
	uniform bool bEnableShadows
	)
{
	vNormalCS = normalize(vNormalCS);

	if (bEnableShadows)
	{
		float4 vPosSM = mul(g_mShadowMap, float4(vPosCS, 1.0f));

		float shade = 1.0;
		if (vPosSM.w > 0.0f)
		{
			vPosSM.xyz /= vPosSM.w;

			float fShadowDepth = g_txShadowMap.SampleLevel(samLinearBorder, vPosSM.xy, 0.0f);

#ifdef SOFT_SHADOWS
			shade = PCSS_Shadow(vPosSM.xy, vPosSM.z, DepthGradient(vPosSM.xy, vPosSM.z), -mul(g_mInvViewLight, float4(vPosCS, 1.0f)).z); // Note the minus sign in before last argument: PCSS_Shadow expects a positive light-space depth (i.e., LHS)
#else
			shade = float(vPosSM.z < fShadowDepth);
#endif
		}
		float dimming = 0.6;
		float diffuseShade = dimming + (1.0f - dimming) * shade;

		vColorOut = Phong(vPosCS, vNormalCS, g_vLightPosCS, g_vAmbientColor, g_vDiffuseColor * diffuseShade, g_vSpecularColor * shade, 30.0f);
	}
	else
	{
		vColorOut = Phong(vPosCS, vNormalCS, g_vLightPosCS, g_vAmbientColor, g_vDiffuseColor, g_vSpecularColor, 30.0f);
	}
}


void vsRenderMeshWithCutSurfaces(
	float3 vPos : POSITION,
	float3 vNormal : NORMAL,
	int iCutID : CUT_ID,
	out float4 vPosOut : SV_POSITION,
	out float3 vPosCSOut : POS_CS,
	out float3 vNormalCSOut : NORMAL_CS,
	out int iCutIDOut : CUT_ID
	)
{
	vPosOut = mul(g_mWorldView, float4(vPos, 1.0f));
	vPosCSOut = vPosOut.xyz;
	vPosOut = mul(g_mProjection, vPosOut);
	vNormalCSOut = normalize( mul(g_mTransInvWorldView, float4(vNormal, 0.0f)).xyz );
	iCutIDOut = iCutID;
}


void psRenderMeshWithCutSurfaces(
	float4 vPos : SV_POSITION,
	float3 vPosCS : POS_CS,
	float3 vNormalCS : NORMAL_CS,
	nointerpolation int iCutID : CUT_ID,
	out float4 vColorOut : SV_TARGET,
	uniform bool bEnableShadows
	)
{
	vNormalCS = normalize(vNormalCS);

	if (bEnableShadows)
	{
		float4 vPosSM = mul(g_mShadowMap, float4(vPosCS, 1.0f));
		vPosSM.xyz /= vPosSM.w;

		float fShadowDepth = g_txShadowMap.SampleLevel(samLinearBorder, vPosSM.xy, 0.0f);

#ifdef SOFT_SHADOWS
		float shade = PCSS_Shadow(vPosSM.xy, vPosSM.z, DepthGradient(vPosSM.xy, vPosSM.z), -mul(g_mInvViewLight, float4(vPosCS, 1.0f)).z); // Note the minus sign in before last argument: PCSS_Shadow expects a positive light-space depth (i.e., LHS)
#else
		float shade = float(vPosSM.z < fShadowDepth);
#endif
		float dimming = 0.6;
		float diffuseShade = dimming + (1.0f - dimming) * shade;

		if (iCutID == 0)
		{
			vColorOut = Phong(vPosCS, vNormalCS, g_vLightPosCS, g_vAmbientColor, g_vDiffuseColor * diffuseShade, g_vSpecularColor * shade, 30.0f);
		}
		else
		{
			vColorOut = Phong(vPosCS, vNormalCS, g_vLightPosCS, g_vAmbientColorCutSurfaces, g_vDiffuseColorCutSurfaces * diffuseShade, g_vSpecularColorCutSurfaces * shade, 30.0f);
		}
	}
	else
	{
		if (iCutID == 0)
		{
			vColorOut = Phong(vPosCS, vNormalCS, g_vLightPosCS, g_vAmbientColor, g_vDiffuseColor, g_vSpecularColor, 30.0f);
		}
		else
		{
			vColorOut = Phong(vPosCS, vNormalCS, g_vLightPosCS, g_vAmbientColorCutSurfaces, g_vDiffuseColorCutSurfaces, g_vSpecularColorCutSurfaces, 30.0f);
		}
	}
}


void vsRenderMeshWithCutSurfacesAndClipping(
	float3 vPos : POSITION,
	float3 vNormal : NORMAL,
	int iCutID : CUT_ID,
	out float4 vPosOut : SV_POSITION,
	out float3 vPosWSOut : POS_WS,
	out float3 vPosCSOut : POS_CS,
	out float3 vNormalCSOut : NORMAL_CS,
	out int iCutIDOut : CUT_ID
	)
{
	vPosOut = mul(g_mWorldView, float4(vPos, 1.0f));
	vPosCSOut = vPosOut.xyz;
	vPosOut = mul(g_mProjection, vPosOut);
	vPosWSOut = mul(g_mWorld, float4(vPos, 1.0f)).xyz;
	vNormalCSOut = normalize( mul(g_mTransInvWorldView, float4(vNormal, 0.0f)).xyz );
	iCutIDOut = iCutID;
}


void psRenderMeshWithCutSurfacesAndClipping(
	float4 vPos : SV_POSITION,
	float3 vPosWS : POS_WS,
	float3 vPosCS : POS_CS,
	float3 vNormalCS : NORMAL_CS,
	int iCutID : CUT_ID,
	out float4 vColorOut : SV_TARGET,
	uniform bool bEnableShadows
	)
{
	if (vPosWS.y < g_fClippingY)
	{
		discard;
	}

	vNormalCS = normalize(vNormalCS);

	if (bEnableShadows)
	{
		float4 vPosSM = mul(g_mShadowMap, float4(vPosCS, 1.0f));
		vPosSM.xyz /= vPosSM.w;

		float fShadowDepth = g_txShadowMap.SampleLevel(samLinearBorder, vPosSM.xy, 0.0f);

#ifdef SOFT_SHADOWS
		float shade = PCSS_Shadow(vPosSM.xy, vPosSM.z, DepthGradient(vPosSM.xy, vPosSM.z), -mul(g_mInvViewLight, float4(vPosCS, 1.0f)).z); // Note the minus sign in before last argument: PCSS_Shadow expects a positive light-space depth (i.e., LHS)
#else
		float shade = float(vPosSM.z < fShadowDepth);
#endif
		float dimming = 0.6;
		float diffuseShade = dimming + (1.0f - dimming) * shade;

		//float fFrontFace = (vNormalCS.z >= 0.0f ? 1.0f : 0.0f);
		if (iCutID == 0)
		{
			vColorOut = Phong(vPosCS, vNormalCS, g_vLightPosCS, g_vAmbientColor, g_vDiffuseColor * diffuseShade, g_vSpecularColor * shade, 30.0f);
		}
		else
		{
			vColorOut = Phong(vPosCS, vNormalCS, g_vLightPosCS, g_vAmbientColorCutSurfaces, g_vDiffuseColorCutSurfaces * diffuseShade, g_vSpecularColorCutSurfaces * shade, 30.0f);
		}
	}
	else
	{
		vColorOut = Phong(vPosCS, vNormalCS, g_vLightPosCS, g_vAmbientColor, g_vDiffuseColor, g_vSpecularColor, 30.0f);
	}
}


void vsRenderMeshIntoShadowMap(
	float3 vPos : POSITION,
	out float4 vPosOut : SV_POSITION
	)
{
	vPosOut = mul(g_mWorldLightSMProjection, float4(vPos, 1.0f));
}


void vsRenderMeshWithFog(
	float3 vPos : POSITION,
	float3 vNormal : NORMAL,
	out float4 vPosOut : SV_POSITION,
	out float3 vPosWSOut : POS_WS,
	out float3 vPosCSOut : POS_CS,
	out float3 vNormalCSOut : NORMAL_CS
	)
{
	vPosOut = mul(g_mWorldView, float4(vPos, 1.0f));
	vPosCSOut = vPosOut.xyz;
	vPosOut = mul(g_mProjection, vPosOut);
	vPosWSOut = mul(g_mWorld, float4(vPos, 1.0f)).xyz;
	vNormalCSOut = normalize( mul(g_mTransInvWorldView, float4(vNormal, 0.0f)).xyz );
}


void psRenderMeshWithFog(
	float4 vPos : SV_POSITION,
	float3 vPosWS : POS_WS,
	float3 vPosCS : POS_CS,
	float3 vNormalCS : NORMAL_CS,
	out float4 vColorOut : SV_TARGET,
	uniform bool bEnableShadows
	)
{
	vNormalCS = normalize(vNormalCS);

	if (bEnableShadows)
	{
		float shade = 1.0;
		float4 vPosSM = mul(g_mShadowMap, float4(vPosCS, 1.0f));
		if (vPosSM.w > 0.0f)
		{
			vPosSM.xyz /= vPosSM.w;
			vPosSM.z = saturate(vPosSM.z);

			float fShadowDepth = g_txShadowMap.SampleLevel(samLinearBorder, vPosSM.xy, 0.0f);

#ifdef SOFT_SHADOWS
			shade = PCSS_Shadow(vPosSM.xy, vPosSM.z, DepthGradient(vPosSM.xy, vPosSM.z), -mul(g_mInvViewLight, float4(vPosCS, 1.0f)).z); // Note the minus sign in before last argument: PCSS_Shadow expects a positive light-space depth (i.e., LHS)
#else
			shade = float(vPosSM.z < fShadowDepth);
#endif
		}
		float dimming = 0.6;
		float diffuseShade = dimming + (1.0f - dimming) * shade;
		float ambientShade = diffuseShade;

		vColorOut = Phong(vPosCS, vNormalCS, g_vLightPosCS, g_vAmbientColor * ambientShade, g_vDiffuseColor * diffuseShade, g_vSpecularColor * shade, 30.0f);

		//if (0.0f <= vPosSM.x && vPosSM.x <= 1.0f && 0.0f <= vPosSM.y && vPosSM.y <= 1.0f && vPosSM.w > 0.0f)
		//{
		//	vColorOut = float4(s, 0.0f, 0.0f, 1.0f);//Phong(vPosCS, vNormalCS, g_vLightPosCS, g_vAmbientColor, s * g_vDiffuseColor, s * g_vSpecularColor, 30.0f);
		//}
		//else
		//{
		//	vColorOut = float4(0.0f, 0.0f, 0.0f, 1.0f);
		//}
	}
	else
	{
		vColorOut = Phong(vPosCS, vNormalCS, g_vLightPosCS, g_vAmbientColor, g_vDiffuseColor, g_vSpecularColor, 30.0f);
	}

	float f = exp(-pow(abs(g_fFogDensity * length(vPosWS - g_vSceneCenter)), g_fFogExponent));
	vColorOut = float4(f * vColorOut.xyz + (1.0f - f) * g_vFogColor, 1.0f);
}


void vsRenderSimulationVertices(
	float3 vPos : POSITION,
	float3 vCenter : CENTER,
	int iType : TYPE,
	uint uiObjectID : SV_INSTANCEID,
	out float4 vPosOut : SV_POSITION,
	out float3 vPosCSOut : POS_CS,
	out nointerpolation float3 vCenterCSOut : CENTER_CS,
	out nointerpolation int iTypeOut : TYPE,
	out nointerpolation uint uiObjectIDOut : OBJECT_ID
	)
{
	vPosOut = float4( vCenter.xyz + g_fSphereRadius * (2.0f * vPos.xyz - float3(1.0f, 1.0f, 1.0f)) , 1.0f);
	vPosOut = mul(g_mWorldView, vPosOut);
	vPosCSOut = vPosOut.xyz;
	vPosOut = mul(g_mProjection, vPosOut);
	vCenterCSOut = mul(g_mWorldView, float4(vCenter.xyz, 1.0f)).xyz;
	iTypeOut = iType;
	uiObjectIDOut = (1 << 28) | uiObjectID;
}


void psRenderSimulationVertices(
	float4 vPos : SV_POSITION,
	float3 vPosCS : POS_CS,
	nointerpolation float3 vCenterCS : CENTER_CS,
	nointerpolation int iType : TYPE,
	nointerpolation uint uiObjectID : OBJECT_ID,
	out float4 vColorOut : SV_TARGET,
	out float fDepthOut : SV_DEPTH,
	out uint uiObjectIDOut : SV_TARGET1
	)
{
	float a = dot(vPosCS, vPosCS); // Ray direction = vPosCS
	float b = -2.0f * dot(vPosCS, vCenterCS);
	float c = dot(vCenterCS, vCenterCS) - g_fSphereRadius * g_fSphereRadius;

	float D = b * b - 4.0f * a * c;
	if (D < 0.0f) { discard; }
	float q = (b + (b >= 0.0f ? 1.0f : -1.0f) * sqrt(D)) / (-2.0f);
	float t = min(q / a, c / q);
	
	float3 vIntersectionCS = t * vPosCS;
	float3 vNormalCS = normalize(vIntersectionCS - vCenterCS);

	float3 vColor = (iType == 1 ? float3(0.2f, 0.0f, 0.0f) : float3(0.0f, 0.2f, 0.0f));
	
	vColorOut = Phong(vPosCS, vNormalCS, g_vLightPosCS, vColor, vColor, vColor, 30.0f);
	
	fDepthOut = -g_mProjection[2][2] - g_mProjection[2][3] / vIntersectionCS.z;

	uiObjectIDOut = uiObjectID;
}


void vsRenderForceArrows(
	float3 vPos : POSITION,
	float3 vNormal : NORMAL,
	int iVertex : VERTEX,
	row_major float4x4 mWorldView : WORLD_VIEW,
	row_major float4x4 mTransInvWorldView : TRANS_INV_WORLD_VIEW,
	out float4 vPosOut : SV_POSITION,
	out float3 vPosCSOut : POS_CS,
	out float3 vNormalCSOut : NORMAL_CS,
	out nointerpolation uint uiObjectIDOut : OBJECT_ID
	)
{
	vPosOut = mul(mWorldView, float4(vPos, 1.0f));
	vPosCSOut = vPosOut.xyz;
	vPosOut = mul(g_mProjection, vPosOut);
	vNormalCSOut = mul(mTransInvWorldView, float4(vNormal, 0.0f)).xyz;
	uiObjectIDOut = (2 << 28) | iVertex;
}


void psRenderForceArrows(
	float4 vPos : SV_POSITION,
	float3 vPosCS : POS_CS,
	float3 vNormalCS : NORMAL_CS,
	nointerpolation uint uiObjectID : OBJECT_ID,
	out float4 vColorOut : SV_TARGET,
	out uint uiObjectIDOut : SV_TARGET1
	)
{
	vNormalCS = normalize(vNormalCS);
	
	float3 vColor = float3(0.3f, 0.3f, 0.0f);
	
	vColorOut = Phong(vPosCS, vNormalCS, g_vLightPosCS, vColor, vColor, vColor, 30.0f);
	
	uiObjectIDOut = uiObjectID;
}


void vsRenderSurfaceVertices(
	float3 vPos : POSITION,
	float3 vCenter : CENTER,
	int iType : TYPE,
	uint uiObjectID : SV_INSTANCEID,
	out float4 vPosOut : SV_POSITION,
	out float3 vPosCSOut : POS_CS,
	out nointerpolation float3 vCenterCSOut : CENTER_CS,
	out nointerpolation int iTypeOut : TYPE,
	out nointerpolation uint uiObjectIDOut : OBJECT_ID
	)
{
	vPosOut = float4( vCenter.xyz + g_fSphereRadius * (2.0f * vPos.xyz - float3(1.0f, 1.0f, 1.0f)) , 1.0f);
	vPosOut = mul(g_mWorldView, vPosOut);
	vPosCSOut = vPosOut.xyz;
	vPosOut = mul(g_mProjection, vPosOut);
	vCenterCSOut = mul(g_mWorldView, float4(vCenter.xyz, 1.0f)).xyz;
	iTypeOut = iType;
	uiObjectIDOut = (3 << 28) | uiObjectID;
}


void psRenderSurfaceVertices(
	float4 vPos : SV_POSITION,
	float3 vPosCS : POS_CS,
	nointerpolation float3 vCenterCS : CENTER_CS,
	nointerpolation int iType : TYPE,
	nointerpolation uint uiObjectID : OBJECT_ID,
	out float4 vColorOut : SV_TARGET,
	out float fDepthOut : SV_DEPTH,
	out uint uiObjectIDOut : SV_TARGET1
	)
{
	float a = dot(vPosCS, vPosCS); // Ray direction = vPosCS
	float b = -2.0f * dot(vPosCS, vCenterCS);
	float c = dot(vCenterCS, vCenterCS) - g_fSphereRadius * g_fSphereRadius;

	float D = b * b - 4.0f * a * c;
	if (D < 0.0f) { discard; }
	float q = (b + (b >= 0.0f ? 1.0f : -1.0f) * sqrt(D)) / (-2.0f);
	float t = min(q / a, c / q);
	
	float3 vIntersectionCS = t * vPosCS;
	float3 vNormalCS = normalize(vIntersectionCS - vCenterCS);

	float3 vColor = float3(iType == 0 || iType == 4 ? 0.2f : 0.0f, iType == 1 || iType == 4 ? 0.2f : 0.0f, iType == 2 ? 0.2f : 0.0f);
	
	vColorOut = Phong(vPosCS, vNormalCS, g_vLightPosCS, vColor, vColor, vColor, 30.0f);
	
	fDepthOut = -g_mProjection[2][2] - g_mProjection[2][3] / vIntersectionCS.z;

	uiObjectIDOut = uiObjectID;
}


technique10 CutSceneRenderer
{
	pass SetSolidRenderMode
	{
		SetRasterizerState( rsCullBack );
	}

	pass SetSolidRenderModeCullNone
	{
		SetRasterizerState( rsCullNone );
	}

	pass SetWireFrameRenderMode
	{
		SetRasterizerState( rsCullNoneWireFrame );
	}

	pass RenderMeshWithFaceNormals			// Renders a mesh with face normals
	{
		SetVertexShader( CompileShader(vs_4_0, vsRenderMeshWithFaceNormals()) );
		SetGeometryShader( CompileShader(gs_4_0, gsRenderMeshWithFaceNormals()) );
		SetPixelShader( CompileShader(ps_4_0, psRenderMesh(false)) );
		
		SetDepthStencilState( dsEnableDepth, 0 );
		SetBlendState( bsDefault, float4( 0.0f, 0.0f, 0.0f, 0.0f ), 0xFFFFFFFF );
	}

	pass RenderMeshWithFaceNormalsAndShadows			// Renders a mesh with face normals and shadows
	{
		SetVertexShader( CompileShader(vs_4_0, vsRenderMeshWithFaceNormals()) );
		SetGeometryShader( CompileShader(gs_4_0, gsRenderMeshWithFaceNormals()) );
		SetPixelShader( CompileShader(ps_4_0, psRenderMesh(true)) );
		
		SetDepthStencilState( dsEnableDepth, 0 );
		SetBlendState( bsDefault, float4( 0.0f, 0.0f, 0.0f, 0.0f ), 0xFFFFFFFF );
	}

	pass RenderMesh							// Renders a mesh with vertex normals
	{
		SetVertexShader( CompileShader(vs_4_0, vsRenderMesh()) );
		SetGeometryShader( NULL );
		SetPixelShader( CompileShader(ps_4_0, psRenderMesh(false)) );
		
		SetDepthStencilState( dsEnableDepth, 0 );
		SetBlendState( bsDefault, float4( 0.0f, 0.0f, 0.0f, 0.0f ), 0xFFFFFFFF );
	}

	pass RenderMeshWithShadows				// Renders a mesh with vertex normals and shadows
	{
		SetVertexShader( CompileShader(vs_4_0, vsRenderMesh()) );
		SetGeometryShader( NULL );
		SetPixelShader( CompileShader(ps_4_0, psRenderMesh(true)) );
		
		SetDepthStencilState( dsEnableDepth, 0 );
		SetBlendState( bsDefault, float4( 0.0f, 0.0f, 0.0f, 0.0f ), 0xFFFFFFFF );
	}

	pass RenderMeshWithCutSurfaces			// Renders a mesh with vertex normals and cut surfaces
	{
		SetVertexShader( CompileShader(vs_4_0, vsRenderMeshWithCutSurfaces()) );
		SetGeometryShader( NULL );
		SetPixelShader( CompileShader(ps_4_0, psRenderMeshWithCutSurfaces(false)) );
		
		SetDepthStencilState( dsEnableDepth, 0 );
		SetBlendState( bsDefault, float4( 0.0f, 0.0f, 0.0f, 0.0f ), 0xFFFFFFFF );
	}

	pass RenderMeshWithCutSurfacesAndShadows			// Renders a mesh with vertex normals, cut surfaces and shadows
	{
		SetVertexShader( CompileShader(vs_4_0, vsRenderMeshWithCutSurfaces()) );
		SetGeometryShader( NULL );
		SetPixelShader( CompileShader(ps_4_0, psRenderMeshWithCutSurfaces(true)) );
		
		SetDepthStencilState( dsEnableDepth, 0 );
		SetBlendState( bsDefault, float4( 0.0f, 0.0f, 0.0f, 0.0f ), 0xFFFFFFFF );
	}

	pass RenderMeshWithCutSurfacesAndClippingAndShadows	// Renders a mesh with vertex normals, cut surfaces, clipping and shadows
	{
		SetVertexShader( CompileShader(vs_4_0, vsRenderMeshWithCutSurfacesAndClipping()) );
		SetGeometryShader( NULL );
		SetPixelShader( CompileShader(ps_4_0, psRenderMeshWithCutSurfacesAndClipping(true)) );
		
		SetRasterizerState( rsCullNone );
		SetDepthStencilState( dsEnableDepth, 0 );
		SetBlendState( bsDefault, float4( 0.0f, 0.0f, 0.0f, 0.0f ), 0xFFFFFFFF );
	}

	pass RenderMeshIntoShadowMap			// Renders a mesh into a depth buffer (shadow map)
	{
		SetVertexShader( CompileShader(vs_4_0, vsRenderMeshIntoShadowMap()) );
		SetGeometryShader( NULL );
		SetPixelShader( NULL );
		
		SetRasterizerState( rsCullNonePolygonOffset );
		SetDepthStencilState( dsEnableDepth, 0 );
		SetBlendState( bsDefault, float4( 0.0f, 0.0f, 0.0f, 0.0f ), 0xFFFFFFFF );
	}

	pass RenderMeshWithFog					// Renders a mesh with vertex normals and fog
	{
		SetVertexShader( CompileShader(vs_4_0, vsRenderMeshWithFog()) );
		SetGeometryShader( NULL );
		SetPixelShader( CompileShader(ps_4_0, psRenderMeshWithFog(false)) );
		
		SetDepthStencilState( dsEnableDepth, 0 );
		SetBlendState( bsDefault, float4( 0.0f, 0.0f, 0.0f, 0.0f ), 0xFFFFFFFF );
	}

	pass RenderMeshWithFogAndShadows		// Renders a mesh with vertex normals, fog and shadows
	{
		SetVertexShader( CompileShader(vs_4_0, vsRenderMeshWithFog()) );
		SetGeometryShader( NULL );
		SetPixelShader( CompileShader(ps_4_0, psRenderMeshWithFog(true)) );
		
		SetDepthStencilState( dsEnableDepth, 0 );
		SetBlendState( bsDefault, float4( 0.0f, 0.0f, 0.0f, 0.0f ), 0xFFFFFFFF );
	}

	pass RenderSimulationVertices			// Renders spheres using ray-casting
	{
		SetVertexShader( CompileShader(vs_4_0, vsRenderSimulationVertices()) );
		SetGeometryShader( NULL );
		SetPixelShader( CompileShader(ps_4_0, psRenderSimulationVertices()) );
		
		SetRasterizerState( rsCullBack );
		SetDepthStencilState( dsEnableDepth, 0 );
		SetBlendState( bsDefault, float4( 0.0f, 0.0f, 0.0f, 0.0f ), 0xFFFFFFFF );
	}

	pass RenderForceArrows
	{
		SetVertexShader( CompileShader(vs_4_0, vsRenderForceArrows()) );
		SetGeometryShader( NULL );
		SetPixelShader( CompileShader(ps_4_0, psRenderForceArrows()) );
		
		SetRasterizerState( rsCullBack );
		SetDepthStencilState( dsEnableDepth, 0 );
		SetBlendState( bsDefault, float4( 0.0f, 0.0f, 0.0f, 0.0f ), 0xFFFFFFFF );
	}

	pass RenderSurfaceVertices			// Renders spheres using ray-casting
	{
		SetVertexShader( CompileShader(vs_4_0, vsRenderSurfaceVertices()) );
		SetGeometryShader( NULL );
		SetPixelShader( CompileShader(ps_4_0, psRenderSurfaceVertices()) );
		
		SetRasterizerState( rsCullBack );
		SetDepthStencilState( dsEnableDepth, 0 );
		SetBlendState( bsDefault, float4( 0.0f, 0.0f, 0.0f, 0.0f ), 0xFFFFFFFF );
	}

}
