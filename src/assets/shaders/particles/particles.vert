#version 120

attribute vec3 aPos;
attribute vec4 aColor;

uniform float uPointSize;
uniform float uRefDistance;
uniform float uMinSize;
uniform float uMaxSize;
uniform int uUseVertexColor;
uniform vec4 uColor;

varying vec4 vColor;

void main()
{
    vec4 eye = gl_ModelViewMatrix * vec4(aPos, 1.0);
    float dist = max(length(eye.xyz), 1e-4);
    float size = uPointSize * (uRefDistance / dist);
    gl_PointSize = clamp(size, uMinSize, uMaxSize);
    gl_Position = gl_ProjectionMatrix * eye;

    if (uUseVertexColor != 0) {
        vColor = aColor;
    } else {
        vColor = uColor;
    }
}