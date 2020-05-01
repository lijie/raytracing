// 相比ch5, 加一个计算射线球体交点处的法线计算

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

float hit_sphere(vec3 center, float radius, ray r) {
    // 计算一元二次方程
    // b^2 - 4ac
    float a = dot(r.direction, r.direction);
    float b = 2.0 * dot(r.origin - center, r.direction);
    float c = dot((r.origin - center), (r.origin - center)) - radius * radius;
    float discriminant = b * b - 4.0 * a * c;
    if (discriminant < 0.0)
        return -1.0;
    
    return (-(b) - sqrt(discriminant)) / (2.0 * a);
}

// 自顶向下返回一个蓝白混合的颜色
vec3 ray_color(ray r) {
    float t = hit_sphere(vec3(0, 0, -1), 0.5, r);
    if (t > 0.0) {
        // 法向量
        // 交点向量 - 圆心向量 = 圆心至交点处的向量 = 法向量
        vec3 n = normalize(ray_point_at(r, t) - vec3(0, 0, -1));
        // scale to [0, 1]
        return 0.5 * (n + 1.0);
        // return vec3(1.0, 0.0, 0.0);
    }
    // scale t to [0, 1]
    t = 0.5 * (r.direction.y + 1.0);
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
