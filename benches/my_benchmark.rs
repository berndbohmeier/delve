use criterion::{Criterion, black_box, criterion_group, criterion_main};

pub fn criterion_benchmark(c: &mut Criterion) {
    // TODO
    c.bench_function("fib 20", |b| b.iter(|| fibonacci(black_box(20))));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
