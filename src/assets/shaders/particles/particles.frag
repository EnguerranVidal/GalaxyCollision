#version 120

uniform vec4 uColor;

void main()
{
    vec2 c = gl_PointCoord * 2.0 - 1.0;
    float r2 = dot(c, c);
    float a = uColor.a * clamp(1.0 - r2, 0.0, 1.0);
    if (a < 0.01)
    discard;
    gl_FragColor = vec4(uColor.rgb, a);
}