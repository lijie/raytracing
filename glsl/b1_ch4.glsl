// 射线
struct ray {
    vec3 origin;    // 原点
    vec3 direction; // 方向, 单位向量
};

ray ray_new(vec3 origin, vec3 direction) {
    ray r;
    r.origin = origin;
    r.direction = direction;
    return r;
}

vec3 ray_point_at(ray r, float t) {
    return r.origin + r.direction * t;
}

// 自顶向下返回一个蓝白混合的颜色
vec3 ray_color(ray r) {
    // scale t to [0, 1]
    float t = 0.5 * (r.direction.y + 1.0);
    return mix(vec3(1.0, 1.0, 1.0), vec3(0.5, 0.7, 1.0), t);
}

void main()
{
    vec3 lower_left_corner = vec3(-2.0, -1.0, -1.0);
    vec3 horizontal = vec3(4.0, 0.0, 0.0);
    vec3 vertical = vec3(0.0, 2.0, 0.0);
    vec3 origin = vec3(0.0, 0.0, 0.0);

    vec2 uv = gl_FragCoord.xy / iResolution.xy;

    vec3 direction = normalize(lower_left_corner + uv.x * horizontal + uv.y * vertical);

    ray r = ray_new(origin, direction);
    vec3 color = ray_color(r);

    gl_FragColor = vec4(color, 1);
}
