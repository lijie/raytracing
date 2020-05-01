// 射线
struct Ray {
    vec3 origin;    // 原点
    vec3 direction; // 方向, 单位向量
}

void main()
{
    vec2 uv = gl_FragCoord.xy / iResolution.xy;
    gl_FragColor = vec4(uv, 0.2, 1);
}
