#version 120

attribute vec3 aPos;

uniform float uPointSize;
uniform float uRefDistance;
uniform float uMinSize;
uniform float uMaxSize;

void main()
{
    vec4 eye = gl_ModelViewMatrix * vec4(aPos, 1.0);
    float dist = max(length(eye.xyz), 1e-4);
    float size = uPointSize * (uRefDistance / dist);
    gl_PointSize = clamp(size, uMinSize, uMaxSize);
    gl_Position = gl_ProjectionMatrix * eye;
}